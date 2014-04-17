# -*- coding: utf-8

# 3 Spheres have an initial velociries: 0, +5, -5
# Their positions and velocities are checking during free fall
# Checks the correctness of NewtonIntegrator and GravityEngine

## Omega
o=Omega() 

## PhysicalParameters 
Density=2400
frictionAngle=radians(35)
sphereRadius=0.05
tc = 0.001
en = 0.3
et = 0.3


sphereMat=O.materials.append(ViscElMat(density=Density,frictionAngle=frictionAngle,tc=tc,en=en,et=et))


v_down = -5.0
v_up = 5.0
g = -9.81
tolerance = 1e-3

id_0=o.bodies.append(utils.sphere((0.0,0,0),0.2,material=sphereMat))    #The body has no initial vertical Velocity
id_down=o.bodies.append(utils.sphere((1.0,0,0),0.2,material=sphereMat)) #The body has an initial vertical Velocity -5
id_up=o.bodies.append(utils.sphere((2.0,0,0),0.2,material=sphereMat))   #The body has an initial vertical Velocity +5

O.bodies[id_down].state.vel[1] = v_down
O.bodies[id_up].state.vel[1] = v_up

## Engines 
o.engines=[
  ForceResetter(),
  InsertionSortCollider([
    Bo1_Sphere_Aabb(),
    Bo1_Facet_Aabb(),
  ]),
  InteractionLoop(
    [Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom()],
    [Ip2_ViscElMat_ViscElMat_ViscElPhys()],
    [Law2_ScGeom_ViscElPhys_Basic()],
  ),
  NewtonIntegrator(damping=0,gravity=[0,g,0]),
  PyRunner(command='checkPos()',iterPeriod=10000),
]

def checkPos():
  if ((O.bodies[id_0].state.pos[1] - getCurrentPos(0))/O.bodies[id_0].state.pos[1] > tolerance): 
    warningMessagePos (0, O.bodies[id_0].state.pos[1], getCurrentPos(0))
  if ((O.bodies[id_down].state.pos[1] - getCurrentPos(v_down))/O.bodies[id_down].state.pos[1] > tolerance): 
    warningMessagePos (v_down, O.bodies[id_down].state.pos[1], getCurrentPos(0))
  if ((O.bodies[id_up].state.pos[1] - getCurrentPos(v_up))/O.bodies[id_up].state.pos[1] > tolerance): 
    warningMessagePos (v_up, O.bodies[id_up].state.pos[1], getCurrentPos(0))
  if ((O.bodies[id_0].state.vel[1] - getCurrentVel(0))/O.bodies[id_0].state.vel[1] > tolerance): 
    warningMessageVel (0, O.bodies[id_0].state.vel[1], getCurrentPos(0))
  if ((O.bodies[id_down].state.vel[1] - getCurrentVel(v_down))/O.bodies[id_down].state.vel[1] > tolerance): 
    warningMessageVel (v_down, O.bodies[id_down].state.vel[1], getCurrentPos(0))
  if ((O.bodies[id_up].state.vel[1] - getCurrentVel(v_up))/O.bodies[id_up].state.vel[1] > tolerance): 
    warningMessageVel (v_up, O.bodies[id_up].state.vel[1], getCurrentPos(0))
  
def getCurrentPos(inVel=0):
  t = O.time+O.dt
  return inVel*t + g*t*t/2

def getCurrentVel(inVel=0):
  t = O.time+O.dt
  return inVel + g*t

def warningMessagePos(inVel, y_pos, y_pos_need):
  print "The body with the initial velocity %.3f, has an y-position %.3f, but it should be %.3f" % (inVel, y_pos, y_pos_need)
  global resultStatus
  resultStatus+=1
  
def warningMessageVel(inVel, y_vel, y_pos_vel):
  print "The body with the initial velocity %.3f, has an y-velocity %.3f, but it should be %.3f" % (inVel, y_vel, y_pos_vel)
  global resultStatus
  resultStatus+=1

O.dt=0.02*tc
O.saveTmp('init');
O.run(1000000)
O.wait()
