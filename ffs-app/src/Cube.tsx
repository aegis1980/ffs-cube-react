function IX(x:number, y:number, z:number, n: number) {
    return x + y * n + z * n * n 
}

class FfsCube {
    size: number;
    dt : number;
    diff : number;
    visc : number;

    s : number[];
    density : number[];

    vx : number[];
    vy : number[];
    vz : number[];

    vx0 : number[];
    vy0 : number[];
    vz0 : number[];


    constructor(size: number, diffusion : number, viscosity : number, dt: number){
        this.size = size;
        this.dt = dt;
        this.diff = diffusion;
        this.visc = viscosity;

        let size3 : number = size**3;
        this.s = Array(size3);
        this.density = Array(size3);

        this.vx  = Array(size3);
        this.vy  = Array(size3);
        this.vz  = Array(size3);

        this.vx0  = Array(size3);
        this.vy0  = Array(size3);
        this.vz0  = Array(size3);
    }

    static setBnd(b : number, x : number[], N : number)
        {
            for(let j = 1; j < N - 1; j++) {
                for(let i = 1; i < N - 1; i++) {
                    x[IX(i, j, 0 ,N  )] = b === 3 ? -x[IX(i, j, 1 ,N )] : x[IX(i, j, 1  ,N)];
                    x[IX(i, j, N-1,N)] = b === 3 ? -x[IX(i, j, N-2,N)] : x[IX(i, j, N-2,N)];
                }
            }
            for(let k = 1; k < N - 1; k++) {
                for(let i = 1; i < N - 1; i++) {
                    x[IX(i, 0  , k,N)] = b === 2 ? -x[IX(i, 1  , k,N)] : x[IX(i, 1  , k,N)];
                    x[IX(i, N-1, k,N)] = b === 2 ? -x[IX(i, N-2, k,N)] : x[IX(i, N-2, k,N)];
                }
            }
            for(let k = 1; k < N - 1; k++) {
                for(let j = 1; j < N - 1; j++) {
                    x[IX(0  , j, k,N)] = b === 1 ? -x[IX(1  , j, k,N)] : x[IX(1  , j, k,N)];
                    x[IX(N-1, j, k,N)] = b === 1 ? -x[IX(N-2, j, k,N)] : x[IX(N-2, j, k,N)];
                }
            }
            
            x[IX(0, 0, 0,N)]       = 0.33 * (x[IX(1, 0, 0,N)]
                                        + x[IX(0, 1, 0,N)]
                                        + x[IX(0, 0, 1,N)]);
            x[IX(0, N-1, 0,N)]     = 0.33 * (x[IX(1, N-1, 0,N)]
                                        + x[IX(0, N-2, 0,N)]
                                        + x[IX(0, N-1, 1,N)]);
            x[IX(0, 0, N-1,N)]     = 0.33 * (x[IX(1, 0, N-1,N)]
                                        + x[IX(0, 1, N-1,N)]
                                        + x[IX(0, 0, N,N)]);
            x[IX(0, N-1, N-1,N)]   = 0.33 * (x[IX(1, N-1, N-1,N)]
                                        + x[IX(0, N-2, N-1,N)]
                                        + x[IX(0, N-1, N-2,N)]);
            x[IX(N-1, 0, 0,N)]     = 0.33 * (x[IX(N-2, 0, 0,N)]
                                        + x[IX(N-1, 1, 0,N)]
                                        + x[IX(N-1, 0, 1,N)]);
            x[IX(N-1, N-1, 0,N)]   = 0.33 * (x[IX(N-2, N-1, 0,N)]
                                        + x[IX(N-1, N-2, 0,N)]
                                        + x[IX(N-1, N-1, 1,N)]);
            x[IX(N-1, 0, N-1,N)]   = 0.33 * (x[IX(N-2, 0, N-1,N)]
                                        + x[IX(N-1, 1, N-1,N)]
                                        + x[IX(N-1, 0, N-2,N)]);
            x[IX(N-1, N-1, N-1,N)] = 0.33 * (x[IX(N-2, N-1, N-1,N)]
                                        + x[IX(N-1, N-2, N-1,N)]
                                        + x[IX(N-1, N-1, N-2,N)]);
        }

    static linSolve(b : number, x : number[], x0 : number[], a : number, c : number, iter : number, N : number)
    {
        let cRecip : number = 1.0 / c;
        for (let k = 0; k < iter; k++) {
            for (let m = 1; m < N - 1; m++) {
                for (let j = 1; j < N - 1; j++) {
                    for (let i = 1; i < N - 1; i++) {
                        x[IX(i, j, m,N)] =
                            (x0[IX(i, j, m,N)]
                                + a*(    x[IX(i+1, j  , m  ,N)]
                                        +x[IX(i-1, j  , m  ,N)]
                                        +x[IX(i  , j+1, m  ,N)]
                                        +x[IX(i  , j-1, m  ,N)]
                                        +x[IX(i  , j  , m+1,N)]
                                        +x[IX(i  , j  , m-1,N)]
                            )) * cRecip;
                    }
                }
            }
            FfsCube.setBnd(b, x, N);
        }
    }

    static diffuse (b : number, x : number[], x0 : number[], diff : number, dt : number, iter : number, N : number)
    {
        let a : number = dt * diff * (N - 2) * (N - 2);
        FfsCube.linSolve(b, x, x0, a, 1 + 6 * a, iter, N);
    }

    static step(cube : FfsCube){
        let N = cube.size;
        let visc = cube.visc;
        let diff = cube.diff;
        let dt = cube.dt;
        let Vx = cube.vx;
        let Vy = cube.vy;
        let Vz = cube.vz;
        let Vx0 = cube.vx0;
        let Vy0 = cube.vy0;
        let Vz0 = cube.vz0;
        let s = cube.s;
        let density = cube.density;

        FfsCube.diffuse(1, Vx0, Vx, visc, dt, 4, N);
        FfsCube.diffuse(2, Vy0, Vy, visc, dt, 4, N);
        FfsCube.diffuse(3, Vz0, Vz, visc, dt, 4, N);
        
        project(Vx0, Vy0, Vz0, Vx, Vy, 4, N);
        
        advect(1, Vx, Vx0, Vx0, Vy0, Vz0, dt, N);
        advect(2, Vy, Vy0, Vx0, Vy0, Vz0, dt, N);
        advect(3, Vz, Vz0, Vx0, Vy0, Vz0, dt, N);
        
        project(Vx, Vy, Vz, Vx0, Vy0, 4, N);
        
        FfsCube.diffuse(0, s, density, diff, dt, 4, N);
        advect(0, density, s, Vx, Vy, Vz, dt, N);

    }
}


export default FfsCube;