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
  const model : any = useRef()
  const fssCube = new FfsCube(prop.size,0.001,0.001,.0001);
  const maxV = 10;
  fssCube.randomValues(10,10);
  const up_axis = new THREE.Vector3(0, 1, 0);
  const dummy = useMemo(() => new THREE.Object3D(), [])
  const tellTales = useMemo(() => {
    const temp :any[] = []
    const w :number = (prop.size*prop.griddim)/2; //so centred on origin
    for (let i = 0; i < prop.size; i++) {
      for (let j = 0; j < prop.size; j++) {
        for (let k = 0; k < prop.size; k++) {
          const index =  FfsCube.ix(i,j,k,prop.size);
          const x :number = (prop.griddim * i)-w;
          const y :number= (prop.griddim * j)-w;
          const z :number= (prop.griddim * k)-w;
          const vx :number= fssCube.vx[index]
          const vy:number = fssCube.vy[index]
          const vz :number=  fssCube.vz[index]
          const density : number = fssCube.density[index]; 
          temp.push({ x, y, z , vx, vy, vz, density })
        }
      }
    }
    return temp
  }, [])

  // Render-loop
  useFrame((state) => {

    const time = state.clock.getElapsedTime();

    const excitation  = Math.sin(time*2)
    for (let i = 0; i < prop.size; i++){
      for (let j = 0; j < prop.size; j++){
        fssCube.setVelocity(i,j,5,3 * excitation,0,0)
      }
    }
    FfsCube.step(fssCube)

    tellTales.forEach((tellTale, i) => {

      let { x, y, z } = tellTale

      const v = new THREE.Vector3(fssCube.vx[i],fssCube.vy[i],fssCube.vz[i])
      const vl = v.clone().length();
      const vn = v.clone().normalize();
      

      // align axis with v
      dummy.scale.set(1,(prop.griddim*0.8)*(vl/maxV),1)
      dummy.quaternion.setFromUnitVectors(up_axis, vn);

      dummy.position.set(
        x, y, z
      )

      dummy.updateMatrix()
      model.current.setMatrixAt(i, dummy.matrix)
    })
    model.current.instanceMatrix.needsUpdate = true
  })

  return (
    <>
      <instancedMesh ref={model} args={[null, null, prop.size**3]}>
        <cylinderBufferGeometry attach="geometry" args={[0, prop.griddim/10, prop.griddim*0.8]} />
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
          gl.setClearColor('lightpink')
          gl.toneMapping = THREE.ACESFilmicToneMapping
          gl.outputEncoding = THREE.sRGBEncoding
        }}>
        <ambientLight intensity={1} />
        <pointLight position={[100, 100, 100]} intensity={2.2} />
        <pointLight position={[-100, -100, -100]} intensity={5} color="blue" />
        <Simulation size={30} griddim={5}/>
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
