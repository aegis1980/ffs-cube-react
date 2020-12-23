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

    static advect(b : number, d : number[], d0 : number[],  velocX: number[], velocY: number[], velocZ: number[], dt: number, N: number)
    {
        let i0, i1, j0, j1, k0, k1 : number;
        
        let dtx : number = dt * (N - 2);
        let dty : number = dt * (N - 2);
        let dtz : number = dt * (N - 2);
        
        let s0, s1, t0, t1, u0, u1 : number;
        let tmp1, tmp2, tmp3, x, y, z : number;
        
        let Nfloat :number = N;
        let ifloat, jfloat, kfloat : number;
        let i, j, k : number;
        
        for(k = 1, kfloat = 1; k < N - 1; k++, kfloat++) {
            for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
                for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
                    tmp1 = dtx * velocX[IX(i, j, k, N)];
                    tmp2 = dty * velocY[IX(i, j, k, N)];
                    tmp3 = dtz * velocZ[IX(i, j, k, N)];
                    x    = ifloat - tmp1; 
                    y    = jfloat - tmp2;
                    z    = kfloat - tmp3;
                    
                    if(x < 0.5) x = 0.5; 
                    if(x > Nfloat + 0.5) x = Nfloat + 0.5; 
                    i0 = Math.floor(x); 
                    i1 = i0 + 1.0;
                    if(y < 0.5) y = 0.5; 
                    if(y > Nfloat + 0.5) y = Nfloat + 0.5; 
                    j0 = Math.floor(y);
                    j1 = j0 + 1.0; 
                    if(z < 0.5) z = 0.5;
                    if(z > Nfloat + 0.5) z = Nfloat + 0.5;
                    k0 = Math.floor(z);
                    k1 = k0 + 1.0;
                    
                    s1 = x - i0; 
                    s0 = 1.0 - s1; 
                    t1 = y - j0; 
                    t0 = 1.0 - t1;
                    u1 = z - k0;
                    u0 = 1.0 - u1;
                    
                    let i0i :number= i0;
                    let i1i :number= i1;
                    let j0i :number = j0;
                    let j1i :number= j1;
                    let k0i :number= k0;
                    let k1i :number = k1;
                    
                    d[IX(i, j, k, N)] = 
                    
                        s0 * ( t0 * (u0 * d0[IX(i0i, j0i, k0i,N)]
                                    +u1 * d0[IX(i0i, j0i, k1i,N)])
                            +( t1 * (u0 * d0[IX(i0i, j1i, k0i,N)]
                                    +u1 * d0[IX(i0i, j1i, k1i,N)])))
                    +s1 * ( t0 * (u0 * d0[IX(i1i, j0i, k0i,N)]
                                    +u1 * d0[IX(i1i, j0i, k1i,N)])
                            +( t1 * (u0 * d0[IX(i1i, j1i, k0i,N)]
                                    +u1 * d0[IX(i1i, j1i, k1i,N)])));
                }
            }
        }
        FfsCube.setBnd(b, d, N);
    }

    
    static project(
            velocX: number[], 
            velocY: number[], 
            velocZ: number[], 
            p: number[], 
            div: number[], 
            iter:number, 
            N:number
        )
    {
        for (let k = 1; k < N - 1; k++) {
            for (let j = 1; j < N - 1; j++) {
                for (let i = 1; i < N - 1; i++) {
                    div[IX(i, j, k, N)] = -0.5*(
                            velocX[IX(i+1, j  , k ,N )]
                            -velocX[IX(i-1, j  , k ,N )]
                            +velocY[IX(i  , j+1, k ,N )]
                            -velocY[IX(i  , j-1, k ,N )]
                            +velocZ[IX(i  , j  , k+1,N)]
                            -velocZ[IX(i  , j  , k-1,N)]
                        )/N;
                    p[IX(i, j, k,N)] = 0;
                }
            }
        }
        FfsCube.setBnd(0, div, N); 
        FfsCube.setBnd(0, p, N);
        FfsCube.linSolve(0, p, div, 1, 6, iter, N);
        
        for (let k = 1; k < N - 1; k++) {
            for (let j = 1; j < N - 1; j++) {
                for (let i = 1; i < N - 1; i++) {
                    velocX[IX(i, j, k,N)] -= 0.5 * (  p[IX(i+1, j, k,N)]
                                                    -p[IX(i-1, j, k,N)]) * N;
                    velocY[IX(i, j, k,N)] -= 0.5 * (  p[IX(i, j+1, k,N)]
                                                    -p[IX(i, j-1, k,N)]) * N;
                    velocZ[IX(i, j, k,N)] -= 0.5 * (  p[IX(i, j, k+1,N)]
                                                    -p[IX(i, j, k-1,N)]) * N;
                }
            }
        }
        FfsCube.setBnd(1, velocX, N);
        FfsCube.setBnd(2, velocY, N);
        FfsCube.setBnd(3, velocZ, N);
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
        
        FfsCube.project(Vx0, Vy0, Vz0, Vx, Vy, 4, N);
        
        FfsCube.advect(1, Vx, Vx0, Vx0, Vy0, Vz0, dt, N);
        FfsCube.advect(2, Vy, Vy0, Vx0, Vy0, Vz0, dt, N);
        FfsCube.advect(3, Vz, Vz0, Vx0, Vy0, Vz0, dt, N);
        
        FfsCube.project(Vx, Vy, Vz, Vx0, Vy0, 4, N);
        
        FfsCube.diffuse(0, s, density, diff, dt, 4, N);
        FfsCube.advect(0, density, s, Vx, Vy, Vz, dt, N);

    }

    
    addDensity(x: number, y: number, z: number, amount: number)
    {
        self.density[IX(x, y, z, self.size)] += amount;
    }

    addVelocity(FluidCube *cube, int x, int y, int z, float amountX, float amountY, float amountZ)
    {
        int N = cube->size;
        int index = IX(x, y, z);
        
        cube->Vx[index] += amountX;
        cube->Vy[index] += amountY;
        cube->Vz[index] += amountZ;
    }
}


export default FfsCube;