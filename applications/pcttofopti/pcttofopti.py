#!/usr/bin/env python
import os
import argparse
import json
from datetime import datetime

import itk
from itk import PCT as pct
from itk import RTK as rtk

from tqdm import tqdm
import torch

m0 = 938.27208943  # MeV/c2
c = 299.792458  # mm/ns

LMComponentType = itk.ctype('float')
LMPixelType = itk.Vector[LMComponentType, 3]
LMImageType = itk.Image[LMPixelType, 2]

RSPPixelType = itk.ctype('double')
RSPImageType = itk.Image[RSPPixelType, 3]


def print_verbose(verbosity, msg):
    if verbosity:
        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}")


class TOFProtonDataset(torch.utils.data.Dataset):

    def __init__(self, lm, angles, polydeg, step, hull):

        self.pos_in = lm[:, 0, :]
        self.pos_out = lm[:, 1, :]
        self.dir_in = lm[:, 2, :]
        self.dir_out = lm[:, 3, :]

        self.tofs = lm[:, 4, 2]
        self.angles = angles
        self.polydeg = polydeg

        assert len(self.pos_in) == len(self.pos_out) == len(
            self.dir_in) == len(self.dir_out) == len(self.tofs) == len(
                self.angles)
        assert torch.all(self.pos_in[0, 2] == self.pos_in[
            1:, 2])  # all protons have the same upstream w
        assert torch.all(self.pos_out[0, 2] == self.pos_out[
            1:, 2])  # all protons have the same downstream w

        w_in = self.pos_in[0][2]
        w_out = self.pos_out[0][2]
        self.ws = torch.arange(w_in, w_out + step, step)

        self.hull = hull

    def get_ws(self):
        return self.ws

    def __len__(self):
        return len(self.tofs)

    def __getitem__(self, idx):

        if self.hull is None:

            mlp = pct.PolynomialMLPFunction.New()
            mlp.SetPolynomialDegree(self.polydeg)
            mlp.Init(self.pos_in[idx].tolist(), self.pos_out[idx].tolist(),
                     self.dir_in[idx].tolist(), self.dir_out[idx].tolist())
            us, vs = mlp.Evaluate(self.ws.tolist())
            us, vs = torch.tensor(us, device=self.ws.device), torch.tensor(
                vs, device=self.ws.device)

        else:

            hull = rtk.QuadricShape.New()
            hull.SetA(self.hull[0])
            hull.SetB(self.hull[1])
            hull.SetC(self.hull[2])
            hull.SetD(self.hull[3])
            hull.SetE(self.hull[4])
            hull.SetF(self.hull[5])
            hull.SetG(self.hull[6])
            hull.SetH(self.hull[7])
            hull.SetI(self.hull[8])
            hull.SetJ(self.hull[9])

            pos_in = self.pos_in[idx]
            pos_out = self.pos_out[idx]
            dir_in = self.dir_in[idx]
            dir_out = self.dir_out[idx]

            intersect_uppp, d_up, _ = hull.IsIntersectedByRay(
                pos_in.tolist(), dir_in.tolist())
            intersect_down, d_down, _ = hull.IsIntersectedByRay(
                pos_out.tolist(), (-dir_out).tolist())
            if not (intersect_uppp and intersect_down):
                us = torch.linspace(pos_in[0],
                                    pos_out[0],
                                    self.ws.shape[0],
                                    device=self.ws.device)
                vs = torch.linspace(pos_in[1],
                                    pos_out[1],
                                    self.ws.shape[0],
                                    device=self.ws.device)
            else:

                int_up = pos_in + dir_in * d_up
                int_down = pos_out - dir_out * d_down

                us = torch.empty_like(self.ws)
                vs = torch.empty_like(self.ws)

                # Before object
                mask_before = self.ws < int_up[2]
                d_before = (self.ws[mask_before] - pos_in[2]) / (int_up[2] -
                                                                 pos_in[2])
                us[mask_before] = pos_in[0] + dir_in[0] * d_up * d_before
                vs[mask_before] = pos_in[1] + dir_in[1] * d_up * d_before

                # After object
                mask_after = self.ws > int_down[2]
                d_after = (self.ws[mask_after] - pos_out[2]) / (int_down[2] -
                                                                pos_out[2])
                us[mask_after] = pos_out[0] - dir_out[0] * d_down * d_after
                vs[mask_after] = pos_out[1] - dir_out[1] * d_down * d_after

                # In object
                mask_in = torch.logical_not(
                    torch.logical_or(mask_before, mask_after))
                mlp = pct.PolynomialMLPFunction.New()
                mlp.SetPolynomialDegree(self.polydeg)
                mlp.Init(int_up.tolist(), int_down.tolist(), dir_in.tolist(),
                         dir_out.tolist())
                us_in, vs_in = mlp.Evaluate(self.ws[mask_in].tolist())
                us[mask_in], vs[mask_in] = torch.tensor(
                    us_in, device=us.device), torch.tensor(vs_in,
                                                           device=vs.device)

        return self.tofs[idx], self.angles[idx], (us, vs, self.ws)


class TOF(torch.autograd.Function):

    @staticmethod
    def forward(ctx, v, d):
        ctx.save_for_backward(v, d)
        return d / v

    @staticmethod
    def backward(ctx, grad_output):
        v, d = ctx.saved_tensors
        return (-d / v**2) * grad_output, None


class Velocity(torch.autograd.Function):

    @staticmethod
    def forward(ctx, e):
        ctx.save_for_backward(e)
        return c * (e / (e + m0)) * torch.sqrt(1 + 2 * (m0 / e))

    @staticmethod
    def backward(ctx, grad_output):
        e, = ctx.saved_tensors
        dv_de = -c * e * torch.sqrt(1 + 2 * m0 / e) / (
            e + m0)**2 + c * torch.sqrt(1 + 2 * m0 / e) / (e + m0) - c * m0 / (
                e * torch.sqrt(1 + 2 * m0 / e) * (e + m0))
        return grad_output * dv_de


class EnergyLoss(torch.autograd.Function):

    @staticmethod
    def forward(ctx, es, rsps, sp_h2o, d):
        ctx.save_for_backward(sp_h2o, d)
        sp = rsps * sp_h2o
        eloss = sp * d
        es2 = es - eloss
        return es2

    @staticmethod
    def backward(ctx, grad_output):
        sp_h2o, d = ctx.saved_tensors
        grad_rsps = -sp_h2o * d * grad_output
        return None, grad_rsps, None, None


class SPWater(torch.autograd.Function):

    @staticmethod
    def forward(ctx, e, p0, p1, p2, p3, p4, p5):
        ctx.save_for_backward(e)
        ctx.params = p0, p1, p2, p3, p4, p5
        return p0 * torch.exp(-p1 * e) + p2 * torch.exp(
            -p3 * e) + p4 * torch.exp(-p5 * e)

    @staticmethod
    def backward(ctx, grad_output):
        e, = ctx.saved_tensors
        p0, p1, p2, p3, p4, p5 = ctx.params
        dsp_water_de = -p0 * p1 * torch.exp(-e * p1) - p2 * p3 * torch.exp(
            -e * p3) - p4 * p5 * torch.exp(-e * p5)
        return grad_output * dsp_water_de, None, None, None, None, None, None


class RSPInterp(torch.autograd.Function):

    @staticmethod
    def forward(ctx, rsp, rsp_spacing, x, y, z):
        upper_bounds = [rsp.shape[i] * rsp_spacing[i] / 2 for i in range(3)]
        lower_bounds = [-upper_bounds[i] for i in range(3)]

        xs = torch.linspace(lower_bounds[0],
                            upper_bounds[0],
                            rsp.shape[0],
                            device=x.device)
        ys = torch.linspace(lower_bounds[1],
                            upper_bounds[1],
                            rsp.shape[1],
                            device=y.device)
        zs = torch.linspace(lower_bounds[2],
                            upper_bounds[2],
                            rsp.shape[2],
                            device=z.device)

        x_idx = torch.searchsorted(xs, x)
        y_idx = torch.searchsorted(ys, y)
        z_idx = torch.searchsorted(zs, z)

        in_mask = ((0 < x_idx) & (x_idx < xs.shape[0])
                   & (0 < y_idx) & (y_idx < ys.shape[0])
                   & (0 < z_idx) & (z_idx < zs.shape[0]))

        oob_value = 0.
        cs = torch.full((len(x), ), oob_value, device=x.device)

        x_in = x[in_mask]
        y_in = y[in_mask]
        z_in = z[in_mask]

        x_idx_in = x_idx[in_mask]
        y_idx_in = y_idx[in_mask]
        z_idx_in = z_idx[in_mask]

        x0, x1 = xs[x_idx_in - 1], xs[x_idx_in]
        y0, y1 = ys[y_idx_in - 1], ys[y_idx_in]
        z0, z1 = zs[z_idx_in - 1], zs[z_idx_in]
        xd = (x_in - x0) / (x1 - x0)
        yd = (y_in - y0) / (y1 - y0)
        zd = (z_in - z0) / (z1 - z0)
        c000 = rsp[x_idx_in - 1, y_idx_in - 1, z_idx_in - 1]
        c001 = rsp[x_idx_in - 1, y_idx_in - 1, z_idx_in]
        c010 = rsp[x_idx_in - 1, y_idx_in, z_idx_in - 1]
        c011 = rsp[x_idx_in - 1, y_idx_in, z_idx_in]
        c100 = rsp[x_idx_in, y_idx_in - 1, z_idx_in - 1]
        c101 = rsp[x_idx_in, y_idx_in - 1, z_idx_in]
        c110 = rsp[x_idx_in, y_idx_in, z_idx_in - 1]
        c111 = rsp[x_idx_in, y_idx_in, z_idx_in]
        c00 = c000 * (1 - xd) + c100 * xd
        c01 = c001 * (1 - xd) + c101 * xd
        c10 = c010 * (1 - xd) + c110 * xd
        c11 = c011 * (1 - xd) + c111 * xd
        c0 = c00 * (1 - yd) + c10 * yd
        c1 = c01 * (1 - yd) + c11 * yd

        cs[in_mask] = c0 * (1 - zd) + c1 * zd

        ctx.rsp_shape = rsp.shape
        ctx.save_for_backward(in_mask, xd, yd, zd, x_idx_in, y_idx_in,
                              z_idx_in)

        return cs

    @staticmethod
    def backward(ctx, grad_output):
        in_mask, xd, yd, zd, x_idx_in, y_idx_in, z_idx_in = ctx.saved_tensors
        rsp_shape = ctx.rsp_shape

        grad_output_in = grad_output[in_mask]
        drsp = torch.zeros(rsp_shape, device=grad_output.device)

        for xx in [False, True]:
            for yy in [False, True]:
                for zz in [False, True]:

                    x = x_idx_in if xx else x_idx_in - 1
                    y = y_idx_in if yy else y_idx_in - 1
                    z = z_idx_in if zz else z_idx_in - 1
                    idx_lin = z + y * rsp_shape[2] + x * rsp_shape[
                        2] * rsp_shape[1]

                    xd_ = xd if xx else 1 - xd
                    yd_ = yd if yy else 1 - yd
                    zd_ = zd if zz else 1 - zd

                    drsp.put_(idx_lin,
                              grad_output_in * xd_ * yd_ * zd_,
                              accumulate=True)

        return drsp, None, None, None, None


def forward_iter(w, xs, ys, zs, es, rsp, rsp_spacing, sp_fit):

    x, y, z = xs[:, w].contiguous(), ys[:, w].contiguous(), zs[:,
                                                               w].contiguous()
    x2, y2, z2 = xs[:, w + 1], ys[:, w + 1], zs[:, w + 1]
    d = torch.sqrt((x - x2)**2 + (y - y2)**2 + (z - z2)**2)

    vel = Velocity.apply(es)
    tof = TOF.apply(vel, d)

    rsps = RSPInterp.apply(rsp, rsp_spacing, x, y, z)
    # Negative RSP is not physical, so the RSPs are clipped to zero
    rsps[rsps < 0.] = 0.

    sp_h2o = SPWater.apply(es, *sp_fit)
    es2 = EnergyLoss.apply(es, rsps, sp_h2o, d)

    # We give some energy (1 keV) to protons that stopped
    # This way they keep going inside the object, although slowly
    # This should give a high TOF that corresponds to a large error w.r.t. data
    es2[es2 <= 0.] = .001

    return tof, es2


def compute_tof(rsp, rsp_spacing, mlps, angles, sp_fit, e0):

    us, vs, ws = mlps
    assert us.shape == vs.shape

    # Apply rotation
    thetas = -torch.deg2rad(angles).unsqueeze(1)
    xs = us * torch.cos(thetas) + ws * torch.sin(thetas)
    ys = vs
    zs = -us * torch.sin(thetas) + ws * torch.cos(thetas)

    number_of_protons = xs.shape[0]
    es = torch.full((number_of_protons, ), e0)
    tofs = torch.full((number_of_protons, ), 0.)

    for w, _ in enumerate(ws[0, :-1]):
        tof, es2 = torch.utils.checkpoint.checkpoint(forward_iter,
                                                     w,
                                                     xs,
                                                     ys,
                                                     zs,
                                                     es,
                                                     rsp,
                                                     rsp_spacing,
                                                     sp_fit,
                                                     use_reentrant=False)
        tofs += tof
        es = es2

    return tofs


def save_img(img, spacing, filename):
    origin = tuple(-(img.shape[i] * spacing[i] / 2) + spacing[i] / 2
                   for i in range(3))

    img_itk = itk.GetImageFromArray(img)
    img_itk.SetSpacing(spacing)
    img_itk.SetOrigin(origin)
    itk.imwrite(img_itk, filename)


def total_variation(img):
    dx = torch.sum(torch.abs(img[1:, :, :] - img[:-1, :, :]))
    dy = torch.sum(torch.abs(img[:, 1:, :] - img[:, :-1, :]))
    dz = torch.sum(torch.abs(img[:, :, 1:] - img[:, :, :-1]))
    return dx + dy + dz


def check_grad_impl():
    torch.set_default_dtype(torch.float64)

    v = torch.rand(10, requires_grad=True) * 1000
    d = torch.rand(10) * 1000
    torch.autograd.gradcheck(TOF.apply, (v, d))

    e = torch.rand(10, requires_grad=True) * 1000
    torch.autograd.gradcheck(Velocity.apply, (e))

    number_of_protons = 10
    es = torch.rand(number_of_protons) * 200
    rsps = torch.rand(number_of_protons, requires_grad=True)
    sp_h2o = torch.rand(number_of_protons) * 10
    d = torch.rand(number_of_protons) * 10
    torch.autograd.gradcheck(EnergyLoss.apply, (es, rsps, sp_h2o, d))

    e = torch.rand(10, requires_grad=True) * 1000
    params = torch.rand(6) * 1000
    torch.autograd.gradcheck(SPWater.apply, (e, *params))

    number_of_protons = 100
    rsp = torch.rand((10, 10, 10), requires_grad=True) * 2
    rsp_spacing = (10, 10, 10)
    ps = torch.rand((number_of_protons, 3)) + -5 * 10
    torch.autograd.gradcheck(RSPInterp.apply,
                             (rsp, rsp_spacing, ps[0], ps[1], ps[2]))

    torch.set_default_dtype(torch.float32)


def pcttofforward(inputs,
                  sp_fit,
                  output_dir=None,
                  rsp_map=None,
                  support_radius=None,
                  initial_energy=200.,
                  step=1.,
                  polydeg=5,
                  hull=None,
                  number_of_iterations=None,
                  learning_rate=.3,
                  tv=0.,
                  batch_size=10000000,
                  num_workers=0,
                  optimizer='SGD',
                  check_grad=False,
                  gpu=False,
                  verbose=False):

    if check_grad:
        print_verbose(verbose, "Checking gradients…")
        check_grad_impl()

    torch.multiprocessing.set_start_method('spawn', force=True)

    device = torch.device(
        "cuda:0" if torch.cuda.is_available() and gpu else "cpu")
    print_verbose(verbose, f"Will use device {device} for computations")

    with open(sp_fit, 'r', encoding='utf-8') as file:
        sp_fit = json.load(file)
    print_verbose(verbose, f"Read SP fit: {sp_fit}")

    print_verbose(verbose, "Attempt at reading checkpoint…")
    checkpoint_path = os.path.join(output_dir, 'checkpoint.pth')
    try:
        checkpoint = torch.load(checkpoint_path, weights_only=False)
        epoch = checkpoint['epoch']
        rsp = checkpoint['rsp']
        rsp_spacing = checkpoint['rsp_spacing']
        optimizer = checkpoint['optimizer']
        loss = checkpoint['loss']
        rsp_shape = rsp.shape
        print_verbose(verbose, "Checkpoint successfully read!")
    except FileNotFoundError as e:
        print_verbose(verbose, "No checkpoint found, starting from scratch.")

        if rsp_map is None:
            rsp_size = (220, 400, 220)
            rsp_shape = (220, 3, 220)
            rsp_spacing = tuple(rsp_size[i] / rsp_shape[i] for i in range(3))
            rsp = torch.full(rsp_shape, 0., device=device)
        else:
            rsp_itk = itk.imread(rsp_map)
            rsp = torch.from_numpy(itk.GetArrayFromImage(rsp_itk)).to(device)
            rsp_shape = rsp.shape
            rsp_spacing = tuple(rsp_itk.GetSpacing())
        print_verbose(
            verbose,
            f"Initialized RSP map with shape {rsp_shape} and spacing {rsp_spacing}"
        )

        epoch = 1

        if optimizer == 'SGD':
            optimizer = torch.optim.SGD([rsp], lr=learning_rate)
        elif optimizer == 'Adagrad':
            optimizer = torch.optim.Adagrad([rsp], lr=learning_rate)
        elif optimizer == 'Adam':
            optimizer = torch.optim.Adam([rsp], lr=learning_rate)
        else:
            raise NotImplementedError from e

        loss = torch.nn.MSELoss(reduction='sum')

    support_mask = None
    if support_radius is not None:
        print_verbose(verbose, "Generating support mask…")
        xx, yy, zz = tuple(
            torch.linspace(-rsp_shape[i] / 2 * rsp_spacing[i], rsp_shape[i] /
                           2 * rsp_spacing[i], rsp_shape[i]) for i in range(3))
        x, _, z = torch.meshgrid(xx, yy, zz)
        d = torch.sqrt(x**2 + z**2)
        support_mask = d > support_radius

    os.makedirs(output_dir, exist_ok=True)

    rsp.requires_grad = True

    number_of_projections = len(inputs)
    angles = torch.linspace(0, 360 - (360 / number_of_projections),
                            number_of_projections)
    print_verbose(verbose, f"Detected {number_of_projections} angles: {angles}")

    print_verbose(verbose, "Reading data…")
    lm = None
    angles_lm = None
    for i, a in zip(inputs, angles):
        lm_itk_i = itk.imread(i)
        lm_i = torch.from_numpy(itk.GetArrayFromImage(lm_itk_i))
        angles_lm_i = torch.full((len(lm_i), ), a)

        lm = lm_i if lm is None else torch.concatenate((lm, lm_i))
        angles_lm = angles_lm_i if angles_lm is None else torch.concatenate(
            (angles_lm, angles_lm_i))
    dataset = TOFProtonDataset(lm, angles_lm, polydeg, step, hull)
    dataloader = torch.utils.data.DataLoader(dataset,
                                             batch_size=batch_size,
                                             num_workers=num_workers,
                                             in_order=False,
                                             persistent_workers=True)

    torch.set_default_device(device)  # from now on everything must run on GPU

    while (number_of_iterations is None) or epoch <= number_of_iterations:
        print_verbose(verbose, f"Starting iteration {epoch}…")

        optimizer.zero_grad()

        loss_tot = 0.

        for tofs_cpu, angle_cpu, mlp_cpu in tqdm(dataloader):
            angle = angle_cpu.to(device)
            mlp = (mlp_cpu[0].to(device), mlp_cpu[1].to(device),
                   mlp_cpu[2].to(device))

            tofs = compute_tof(rsp, rsp_spacing, mlp, angle, sp_fit,
                               initial_energy)
            tofs_m = tofs_cpu.to(device)

            loss_n = loss(tofs, tofs_m) + tv * total_variation(rsp)
            loss_tot += loss_n.item()

            loss_n.backward()

        if support_mask is not None:
            rsp.grad[support_mask] = 0.

        optimizer.step()

        with torch.no_grad():
            rsp[rsp < 0.] = 0.

        save_img(rsp.grad.cpu(), rsp_spacing,
                 os.path.join(output_dir, f'grad_loss_{epoch}.mhd'))
        save_img(rsp.cpu().detach().numpy(), rsp_spacing,
                 os.path.join(output_dir, f'rsp_{epoch}.mhd'))
        with open(os.path.join(output_dir, 'loss.txt'), 'a',
                  encoding='utf-8') as f:
            f.write(str(loss_tot) + os.linesep)

        epoch += 1

        checkpoint = {
            'epoch': epoch,
            'rsp': rsp,
            'rsp_spacing': rsp_spacing,
            'optimizer': optimizer,
            'loss': loss
        }
        torch.save(checkpoint, checkpoint_path)

    print(rsp.sum())


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i',
                        '--inputs',
                        help="Input list-mode data",
                        required=True,
                        nargs='+')
    parser.add_argument('-o', '--output-dir', help="Output directory")
    parser.add_argument('--sp-fit',
                        help="Fit of stopping power of water",
                        required=True)
    parser.add_argument('-r',
                        '--rsp-map',
                        help="Initial relative stopping power map")
    parser.add_argument('-s',
                        '--support-radius',
                        help="Zero-out gradients outside of object support",
                        type=float,
                        required=False)
    parser.add_argument('-e',
                        '--initial-energy',
                        help="Initial energy of the protons in MeV",
                        type=float,
                        default=200.)
    parser.add_argument('--step',
                        help="Step size (in mm along z)",
                        type=float,
                        default=1.)
    parser.add_argument('--polydeg',
                        help="Polynomial degree for MLP estimation",
                        type=int,
                        default=5)
    parser.add_argument('--hull',
                        help="Quadric representing the hull of the object",
                        type=float,
                        nargs=10,
                        required=False)
    parser.add_argument('--optimizer', help="Optimizer", default='Adam')
    parser.add_argument('-n',
                        '--number-of-iterations',
                        help="Number of iterations",
                        type=int,
                        required=False)
    parser.add_argument('-l',
                        '--learning-rate',
                        help="Learning rate",
                        type=float,
                        default=.3)
    parser.add_argument('--tv',
                        help="Total variation regularization",
                        type=float,
                        default=0.)
    parser.add_argument('-b',
                        '--batch-size',
                        help="Batch size",
                        type=int,
                        required=False,
                        default=100_000)
    parser.add_argument('-w',
                        '--num-workers',
                        help="Number of workers to compute MLPs",
                        type=int,
                        required=False,
                        default=0)
    parser.add_argument('--check-grad',
                        help="Check gradient implementation",
                        default=False,
                        action='store_true')
    parser.add_argument('--gpu',
                        help="Use GPU",
                        default=False,
                        action='store_true')
    parser.add_argument('--verbose',
                        '-v',
                        help="Verbose execution",
                        default=False,
                        action='store_true')
    args_info = parser.parse_args()

    pcttofforward(**vars(args_info))


if __name__ == '__main__':
    main()
