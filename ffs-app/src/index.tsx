import * as THREE from 'three'
import ReactDOM from 'react-dom'
import React, { useRef, useMemo } from 'react'
import { Canvas, useFrame , extend, useThree, ReactThreeFiber } from 'react-three-fiber'
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls'
import  FfsCube from './Cube'

import './styles.css'

extend({ OrbitControls })

declare global {
  namespace JSX {
    interface IntrinsicElements {
      orbitControls: ReactThreeFiber.Object3DNode<OrbitControls, typeof OrbitControls>
    }
  }
}
interface SimProp {
  size: number;
  griddim : number;
}

function Simulation(prop : SimProp)  {
  const mesh : any = useRef()
  const dummy = useMemo(() => new THREE.Object3D(), [])
  const fssCube = new FfsCube(prop.size,1,1,1);

  const particles = useMemo(() => {
    const temp :any[] = []
    for (let i = 0; i < prop.size; i++) {
      for (let j = 0; j < prop.size; j++) {
        for (let k = 0; k < prop.size; k++) {
          const index =  FfsCube.getIndex(i,j,k,prop.size);
          const x :number = prop.griddim * i;
          const y :number= prop.griddim * j;
          const z :number= prop.griddim * k;
          const vx :number= fssCube.vx[0]
          const vy:number = fssCube.vy[0]
          const vz :number=  fssCube.vz[0]
          const density : number = fssCube.density[0]; 
          temp.push({ x, y, z , vx, vy, vz, density })
        }
      }
    }
    return temp
  }, [prop.size*prop.size*prop.size])

  useFrame((state) => {
    FfsCube.step(fssCube)

    particles.forEach((particle, i) => {
      let { x, y, z, vx, vy, vz, density } = particle
      dummy.position.set(
        x, y, z
      )
      dummy.updateMatrix()
      mesh.current.setMatrixAt(i, dummy.matrix)
    })
    mesh.current.instanceMatrix.needsUpdate = true
  })

  return (
    <>
      <instancedMesh ref={mesh} args={[null, null, prop.size**3]}>
        <cylinderBufferGeometry attach="geometry" args={[3, 3, prop.griddim*0.8]} />
        <meshPhongMaterial attach="material" color="blue" />
      </instancedMesh>
    </>
  )
}

function App() {

  const Scene = () => {
    const {
      camera,
      gl: { domElement }
    } = useThree()
    return (
      <>
        <orbitControls args={[camera, domElement]} />
      </>
    )
  }

  return (
    <div style={{ width: '100%', height: '100%' }}>
      <Canvas
        gl={{ alpha: false, antialias: false }}
        camera={{ position: [0, 0, 70], near:5, far:20000 }}
        onCreated={({ gl }) => {
          gl.setClearColor('white')
          gl.toneMapping = THREE.ACESFilmicToneMapping
          gl.outputEncoding = THREE.sRGBEncoding
        }}>
        <ambientLight intensity={1} />
        <pointLight position={[100, 100, 100]} intensity={2.2} />
        <pointLight position={[-100, -100, -100]} intensity={5} color="blue" />
        <Simulation size={10} griddim={150}/>
        <Scene/>
      </Canvas>
    </div>
  )
}

ReactDOM.render(<App />, document.getElementById('root'))


// If you want to start measuring performance in your app, pass a function
// to log results (for example: reportWebVitals(console.log))
// or send to an analytics endpoint. Learn more: https://bit.ly/CRA-vitals
//reportWebVitals();
