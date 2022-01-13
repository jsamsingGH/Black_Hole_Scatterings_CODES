import math
import numpy
from amuse.lab import *
from amuse.couple import bridge
from matplotlib import pyplot
from plummer_potential import Plummer_potential

class CodeWithFriction(bridge.GravityCodeInField):
    
    def kick_with_field_code(self, particles, field_code, dt):
        self.LnL = 3.7

        R = particles.position.length()
        vx = particles.vx.mean()
        vy = particles.vy.mean()
        vz = particles.vz.mean()
        rho = field_code.density(R)
        vc = field_code.circular_velocity(R)
        X = 0.34 
        m = particles.mass.sum()
        ax = -4*numpy.pi*self.LnL*constants.G**2 * rho*m*(vx/vc**3)*X
        ay = -4*numpy.pi*self.LnL*constants.G**2 * rho*m*(vy/vc**3)*X
        az = -4*numpy.pi*self.LnL*constants.G**2 * rho*m*(vz/vc**3)*X
        self.update_velocities(particles, dt, ax, ay, az)

    def drift(self, tend): 
        pass

def get_system_state(time, black_holes):
    t = []
    r = []
    for bi in black_holes:
        t.append(time.value_in(units.Myr))
        r.append(bi.position.length().value_in(units.parsec))
    return t, r
#PROV HER AT UDLAES x,y,z OG PLOT!!!

    
def evolve_cluster_in_potential(gravity, black_holes, t_end, dt,
                                channel_to_framework):

    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot_prev = Etot_init

    x = []
    y = []
    t, r = get_system_state(gravity.model_time, black_holes)
    x.append(t)
    y.append(r)
    while gravity.model_time < t_end:

        gravity.evolve_model(gravity.model_time + dt)
        channel_to_framework.copy()
        
        t, r = get_system_state(gravity.model_time, black_holes)
        x.append(t)
        y.append(r)        
#        x.append(gravity.model_time.value_in(units.Myr))
#        y.append(black_holes.position.length().value_in(units.parsec))

        Etot_prev_se = gravity.kinetic_energy + gravity.potential_energy
        Ekin = gravity.kinetic_energy 
        Epot = gravity.potential_energy
        Etot = Ekin + Epot
        print("T=", gravity.model_time, end=' ')
        print("E= ", Etot.in_(units.erg), "Q= ", Ekin/Epot, end=' ')
        print("dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot) 
        Etot_prev = Etot
    return x, y

def integrate_black_holes_in_potential(black_holes, Plummer, t_end, dt, dt_bridge, converter):

    cluster_gravity = ph4(converter)#Mikkola(converter)
    cluster_gravity.parameters.lightspeed = constants.c

    cluster_gravity.particles.add_particles(black_holes)
    channel_from_gravity_to_framework = cluster_gravity.particles.new_channel_to(black_holes)

    friction_code = CodeWithFriction(cluster_gravity, (Plummer,), do_sync=True, verbose=False, radius_is_eps=False, h_smooth_is_eps=False, zero_smoothing=False)
        
    gravity = bridge.Bridge(use_threading=False)
    gravity.add_system(cluster_gravity, (Plummer,) )
    gravity.add_code(friction_code)
    gravity.timestep = dt_bridge

    x, y = evolve_cluster_in_potential(gravity, black_holes, t_end, dt,
                                       channel_from_gravity_to_framework)
    gravity.stop()
    return x, y

def plot_orbit(x, y):
    fig = pyplot.figure(figsize=(8,8))	
    pyplot.xlabel("t [Myr]")
    pyplot.ylabel("Y [pc]")
    pyplot.plot(x, y, lw=1)
    pyplot.scatter(x[0], y[0], lw=1)
    pyplot.show()

def main(t_end, M_cl, R_Pl, dt_out, dt_bridge):
    dt_bridge = min(dt_bridge, dt_out)
    Plummer = Plummer_potential(M_cl, R_Pl)
    converter=nbody_system.nbody_to_si(M_cl, R_Pl)
    #N = 100
    black_holes = new_plummer_gas_model(20, convert_nbody=converter)
    black_holes.mass = 100 | units.MSun
    black_holes.radius = 2 * constants.G * black_holes.mass/constants.c**2

    x, y = integrate_black_holes_in_potential(black_holes, Plummer, t_end, dt_out, dt_bridge, converter)

    plot_orbit(x, y)

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("--dt_out", dest="dt_out", type="float", default = 0.01|units.Myr,
                      help="output_timestep [%default]")
    result.add_option("--dt_bridge", unit=units.yr,
                      dest="dt_bridge", type="float", default = 1000|units.yr,
                      help="bridge timestep [%default]")
    result.add_option("-M", unit=units.MSun,
                      dest="M_cl", type="float",default = 1.0e+6 | units.MSun,
                      help="cluster mass [%default]")
    result.add_option("-R", unit= units.parsec,
                      dest="R_Pl", type="float",default = 1.0 | units.parsec,
                      help="Plummer radius [%default]")
    result.add_option("-t", unit= units.Myr,
                      dest="t_end", type="float", default = 1.0 | units.Myr,
                      help="end time of the simulation [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

