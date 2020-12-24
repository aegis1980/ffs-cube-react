//import { isConstructSignatureDeclaration } from "typescript";
var FfsCube = /** @class */ (function () {
    function FfsCube(size, diffusion, viscosity, dt) {
        var _this = this;
        this.toString = function () {
            var index = FfsCube.getIndex(_this.size / 2, _this.size / 2, _this.size / 2, _this.size);
            return '' + (_this.vx[index]);
        };
        this.size = size;
        this.dt = dt;
        this.diff = diffusion;
        this.visc = viscosity;
        var size3 = Math.pow(size, 3);
        this.s = Array(size3).fill(0);
        this.density = Array(size3).fill(0);
        this.vx = Array(size3).fill(0);
        this.vy = Array(size3).fill(0);
        this.vz = Array(size3).fill(0);
        this.vx0 = Array(size3).fill(0);
        this.vy0 = Array(size3).fill(0);
        this.vz0 = Array(size3).fill(0);
    }
    FfsCube.setBnd = function (b, x, N) {
        for (var j = 1; j < N - 1; j++) {
            for (var i = 1; i < N - 1; i++) {
                x[FfsCube.getIndex(i, j, 0, N)] = b === 3 ? -x[FfsCube.getIndex(i, j, 1, N)] : x[FfsCube.getIndex(i, j, 1, N)];
                x[FfsCube.getIndex(i, j, N - 1, N)] = b === 3 ? -x[FfsCube.getIndex(i, j, N - 2, N)] : x[FfsCube.getIndex(i, j, N - 2, N)];
            }
        }
        for (var k = 1; k < N - 1; k++) {
            for (var i = 1; i < N - 1; i++) {
                x[FfsCube.getIndex(i, 0, k, N)] = b === 2 ? -x[FfsCube.getIndex(i, 1, k, N)] : x[FfsCube.getIndex(i, 1, k, N)];
                x[FfsCube.getIndex(i, N - 1, k, N)] = b === 2 ? -x[FfsCube.getIndex(i, N - 2, k, N)] : x[FfsCube.getIndex(i, N - 2, k, N)];
            }
        }
        for (var k = 1; k < N - 1; k++) {
            for (var j = 1; j < N - 1; j++) {
                x[FfsCube.getIndex(0, j, k, N)] = b === 1 ? -x[FfsCube.getIndex(1, j, k, N)] : x[FfsCube.getIndex(1, j, k, N)];
                x[FfsCube.getIndex(N - 1, j, k, N)] = b === 1 ? -x[FfsCube.getIndex(N - 2, j, k, N)] : x[FfsCube.getIndex(N - 2, j, k, N)];
            }
        }
        x[FfsCube.getIndex(0, 0, 0, N)] = 0.33 * (x[FfsCube.getIndex(1, 0, 0, N)]
            + x[FfsCube.getIndex(0, 1, 0, N)]
            + x[FfsCube.getIndex(0, 0, 1, N)]);
        x[FfsCube.getIndex(0, N - 1, 0, N)] = 0.33 * (x[FfsCube.getIndex(1, N - 1, 0, N)]
            + x[FfsCube.getIndex(0, N - 2, 0, N)]
            + x[FfsCube.getIndex(0, N - 1, 1, N)]);
        x[FfsCube.getIndex(0, 0, N - 1, N)] = 0.33 * (x[FfsCube.getIndex(1, 0, N - 1, N)]
            + x[FfsCube.getIndex(0, 1, N - 1, N)]
            + x[FfsCube.getIndex(0, 0, N, N)]);
        x[FfsCube.getIndex(0, N - 1, N - 1, N)] = 0.33 * (x[FfsCube.getIndex(1, N - 1, N - 1, N)]
            + x[FfsCube.getIndex(0, N - 2, N - 1, N)]
            + x[FfsCube.getIndex(0, N - 1, N - 2, N)]);
        x[FfsCube.getIndex(N - 1, 0, 0, N)] = 0.33 * (x[FfsCube.getIndex(N - 2, 0, 0, N)]
            + x[FfsCube.getIndex(N - 1, 1, 0, N)]
            + x[FfsCube.getIndex(N - 1, 0, 1, N)]);
        x[FfsCube.getIndex(N - 1, N - 1, 0, N)] = 0.33 * (x[FfsCube.getIndex(N - 2, N - 1, 0, N)]
            + x[FfsCube.getIndex(N - 1, N - 2, 0, N)]
            + x[FfsCube.getIndex(N - 1, N - 1, 1, N)]);
        x[FfsCube.getIndex(N - 1, 0, N - 1, N)] = 0.33 * (x[FfsCube.getIndex(N - 2, 0, N - 1, N)]
            + x[FfsCube.getIndex(N - 1, 1, N - 1, N)]
            + x[FfsCube.getIndex(N - 1, 0, N - 2, N)]);
        x[FfsCube.getIndex(N - 1, N - 1, N - 1, N)] = 0.33 * (x[FfsCube.getIndex(N - 2, N - 1, N - 1, N)]
            + x[FfsCube.getIndex(N - 1, N - 2, N - 1, N)]
            + x[FfsCube.getIndex(N - 1, N - 1, N - 2, N)]);
    };
    FfsCube.linSolve = function (b, x, x0, a, c, iter, N) {
        var cRecip = 1.0 / c;
        for (var k = 0; k < iter; k++) {
            for (var m = 1; m < N - 1; m++) {
                for (var j = 1; j < N - 1; j++) {
                    for (var i = 1; i < N - 1; i++) {
                        x[FfsCube.getIndex(i, j, m, N)] =
                            (x0[FfsCube.getIndex(i, j, m, N)]
                                + a * (x[FfsCube.getIndex(i + 1, j, m, N)]
                                    + x[FfsCube.getIndex(i - 1, j, m, N)]
                                    + x[FfsCube.getIndex(i, j + 1, m, N)]
                                    + x[FfsCube.getIndex(i, j - 1, m, N)]
                                    + x[FfsCube.getIndex(i, j, m + 1, N)]
                                    + x[FfsCube.getIndex(i, j, m - 1, N)])) * cRecip;
                    }
                }
            }
            FfsCube.setBnd(b, x, N);
        }
    };
    FfsCube.diffuse = function (b, x, x0, diff, dt, iter, N) {
        var a = dt * diff * (N - 2) * (N - 2);
        FfsCube.linSolve(b, x, x0, a, 1 + 6 * a, iter, N);
    };
    FfsCube.advect = function (b, d, d0, velocX, velocY, velocZ, dt, N) {
        var i0, i1, j0, j1, k0, k1;
        var dtx = dt * (N - 2);
        var dty = dt * (N - 2);
        var dtz = dt * (N - 2);
        var s0, s1, t0, t1, u0, u1;
        var tmp1, tmp2, tmp3, x, y, z;
        var Nfloat = N;
        var ifloat, jfloat, kfloat;
        var i, j, k;
        for (k = 1, kfloat = 1; k < N - 1; k++, kfloat++) {
            for (j = 1, jfloat = 1; j < N - 1; j++, jfloat++) {
                for (i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
                    tmp1 = dtx * velocX[FfsCube.getIndex(i, j, k, N)];
                    tmp2 = dty * velocY[FfsCube.getIndex(i, j, k, N)];
                    tmp3 = dtz * velocZ[FfsCube.getIndex(i, j, k, N)];
                    x = ifloat - tmp1;
                    y = jfloat - tmp2;
                    z = kfloat - tmp3;
                    if (x < 0.5)
                        x = 0.5;
                    if (x > Nfloat + 0.5)
                        x = Nfloat + 0.5;
                    i0 = Math.floor(x);
                    i1 = i0 + 1.0;
                    if (y < 0.5)
                        y = 0.5;
                    if (y > Nfloat + 0.5)
                        y = Nfloat + 0.5;
                    j0 = Math.floor(y);
                    j1 = j0 + 1.0;
                    if (z < 0.5)
                        z = 0.5;
                    if (z > Nfloat + 0.5)
                        z = Nfloat + 0.5;
                    k0 = Math.floor(z);
                    k1 = k0 + 1.0;
                    s1 = x - i0;
                    s0 = 1.0 - s1;
                    t1 = y - j0;
                    t0 = 1.0 - t1;
                    u1 = z - k0;
                    u0 = 1.0 - u1;
                    var i0i = i0;
                    var i1i = i1;
                    var j0i = j0;
                    var j1i = j1;
                    var k0i = k0;
                    var k1i = k1;
                    d[FfsCube.getIndex(i, j, k, N)] =
                        s0 * (t0 * (u0 * d0[FfsCube.getIndex(i0i, j0i, k0i, N)]
                            + u1 * d0[FfsCube.getIndex(i0i, j0i, k1i, N)])
                            + (t1 * (u0 * d0[FfsCube.getIndex(i0i, j1i, k0i, N)]
                                + u1 * d0[FfsCube.getIndex(i0i, j1i, k1i, N)])))
                            + s1 * (t0 * (u0 * d0[FfsCube.getIndex(i1i, j0i, k0i, N)]
                                + u1 * d0[FfsCube.getIndex(i1i, j0i, k1i, N)])
                                + (t1 * (u0 * d0[FfsCube.getIndex(i1i, j1i, k0i, N)]
                                    + u1 * d0[FfsCube.getIndex(i1i, j1i, k1i, N)])));
                }
            }
        }
        FfsCube.setBnd(b, d, N);
    };
    FfsCube.project = function (velocX, velocY, velocZ, p, div, iter, N) {
        for (var k = 1; k < N - 1; k++) {
            for (var j = 1; j < N - 1; j++) {
                for (var i = 1; i < N - 1; i++) {
                    div[FfsCube.getIndex(i, j, k, N)] = -0.5 * (velocX[FfsCube.getIndex(i + 1, j, k, N)]
                        - velocX[FfsCube.getIndex(i - 1, j, k, N)]
                        + velocY[FfsCube.getIndex(i, j + 1, k, N)]
                        - velocY[FfsCube.getIndex(i, j - 1, k, N)]
                        + velocZ[FfsCube.getIndex(i, j, k + 1, N)]
                        - velocZ[FfsCube.getIndex(i, j, k - 1, N)]) / N;
                    p[FfsCube.getIndex(i, j, k, N)] = 0;
                }
            }
        }
        FfsCube.setBnd(0, div, N);
        FfsCube.setBnd(0, p, N);
        FfsCube.linSolve(0, p, div, 1, 6, iter, N);
        for (var k = 1; k < N - 1; k++) {
            for (var j = 1; j < N - 1; j++) {
                for (var i = 1; i < N - 1; i++) {
                    velocX[FfsCube.getIndex(i, j, k, N)] -= 0.5 * (p[FfsCube.getIndex(i + 1, j, k, N)]
                        - p[FfsCube.getIndex(i - 1, j, k, N)]) * N;
                    velocY[FfsCube.getIndex(i, j, k, N)] -= 0.5 * (p[FfsCube.getIndex(i, j + 1, k, N)]
                        - p[FfsCube.getIndex(i, j - 1, k, N)]) * N;
                    velocZ[FfsCube.getIndex(i, j, k, N)] -= 0.5 * (p[FfsCube.getIndex(i, j, k + 1, N)]
                        - p[FfsCube.getIndex(i, j, k - 1, N)]) * N;
                }
            }
        }
        FfsCube.setBnd(1, velocX, N);
        FfsCube.setBnd(2, velocY, N);
        FfsCube.setBnd(3, velocZ, N);
    };
    FfsCube.getIndex = function (x, y, z, n) {
        return x + y * n + z * n * n;
    };
    FfsCube.step = function (cube) {
        var N = cube.size;
        var visc = cube.visc;
        var diff = cube.diff;
        var dt = cube.dt;
        var Vx = cube.vx;
        var Vy = cube.vy;
        var Vz = cube.vz;
        var Vx0 = cube.vx0;
        var Vy0 = cube.vy0;
        var Vz0 = cube.vz0;
        var s = cube.s;
        var density = cube.density;
        FfsCube.diffuse(1, Vx0, Vx, visc, dt, 4, N);
        FfsCube.diffuse(2, Vy0, Vy, visc, dt, 4, N);
        FfsCube.diffuse(3, Vz0, Vz, visc, dt, 4, N);
        FfsCube.project(Vx0, Vy0, Vz0, Vx, Vy, 4, N);
        FfsCube.advect(1, Vx, Vx0, Vx0, Vy0, Vz0, dt, N);
        FfsCube.advect(2, Vy, Vy0, Vx0, Vy0, Vz0, dt, N);
        FfsCube.advect(3, Vz, Vz0, Vx0, Vy0, Vz0, dt, N);
        FfsCube.project(Vx, Vy, Vz, Vx0, Vy0, 4, N);
        FfsCube.diffuse(0, s, density, diff, dt, 4, N);
        FfsCube.advect(0, density, s, Vx, Vy, Vz, dt, N);
    };
    FfsCube.prototype.addDensity = function (x, y, z, amount) {
        this.density[FfsCube.getIndex(x, y, z, this.size)] += amount;
    };
    FfsCube.prototype.addVelocity = function (x, y, z, amountX, amountY, amountZ) {
        var index = FfsCube.getIndex(x, y, z, this.size);
        this.vx[index] += amountX;
        this.vy[index] += amountY;
        this.vz[index] += amountZ;
    };
    return FfsCube;
}());
//if (typeof require !== 'undefined' && require.main === module) {
console.log('Running...');
var dt = 0.1;
var mCube = new FfsCube(10, 0.4, 0.4, dt);
mCube.addDensity(5, 5, 5, 10);
mCube.addVelocity(2, 2, 2, 10, 10, 10);
console.log(mCube.toString());
var t = 0;
while (t < 100) {
    FfsCube.step(mCube);
    t += dt;
}
console.log(mCube.toString());
//}
//export default FfsCube;
