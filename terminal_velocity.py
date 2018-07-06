#!/usr/bin/env python
#
# This Python code calculates the terminal velocity of an object falling
# through a fluid.
#
# Copyright (C) 2018, Huidae Cho <https://geni.isnew.info>
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
################################################################################

# Rain drop problem by Huidae Cho
#
# There are two rain droplets of the same diameter of 4mm falling from 1km and
# 2km, individually. Write Python code to find their terminal velocity when
# they hit the ground and see if the terminal velocity is a function of the
# height.
#
# Assume that there is no wind, the densities of air and rain droplets are 1.2
# kg/m^3 and 998 kg/m^3, respectively, the kinematic viscosity of air is
# 1.515*10^-5 m^2/s, and the rain droplets are a sphere and do not experience
# deformation while descending. Use Stokes' law (Cd = 24/Re for Re < 0.5) and
# Clift and Gauvin drag correlation (Cd = 24/Re * (1+0.15*Re^0.687) + 0.42 /
# (1+4.25*10^4*Re^-1.16) for Re < 3*10^5) to calculate the drag coefficient
# where Cd and Re represent the drag coefficient and Reynolds number,
# respectively.

rain_dia = 4e-3 # m
rain_dens = 998 # kg/m^3
air_dens = 1.2 # kg/m^3
air_kine_vis = 1.515e-5 # m^2/s

obj_dia = rain_dia
obj_dens = rain_dens
flu_dens = air_dens
flu_kine_vis = air_kine_vis

# known problem for validation; vel_term should be about 0.438 m/s
#obj_dia = 0.02
#obj_dens = 1.3*998
#flu_kine_vis = 1e-6
#flu_dens = 998

# parameters for method 1
vel_guess = 1 # m/s
vel_thresh = 1e-3 # m/s
max_iter = 1000

# parameters for method 2
time_res = 1e-5 # s
force_bal_thresh = 1e-10 # N

# constants
grav = 9.81 # m/s^2

from math import pi, sqrt

def calc_drag_coef(vel, dia, vis):
    # not a function of height
    rey = vel*dia/vis
    if rey < 0.5:
        drag_coef = 24/rey
    elif rey < 3e5:
        drag_coef = 24/rey*(1+0.15*rey**0.687)+0.42/(1+4.25*1e4*rey**-1.16)
    else:
        raise Exception('Reynolds number too large')
    return drag_coef

def calc_drag_force(drag_coef, area, dens, vel):
    # not a function of height
    return drag_coef*area*dens*vel**2/2

def calc_obj_vel(obj_dia, obj_dens, flu_dens, drag_coef):
    # not a function of height
    return sqrt((obj_dens-flu_dens)*grav*4/3*obj_dia/(drag_coef*flu_dens))

def calc_obj_term_vel_1(vel_guess, vel_thresh, max_iter, obj_dia, obj_dens,
        flu_dens, flu_kine_vis, verbose):
    # not a function of height
    obj_vol = pi*obj_dia**3/6
    obj_area = pi*obj_dia**2/4

    buoy_force = flu_dens*grav*obj_vol
    obj_weight = obj_dens*grav*obj_vol

    found = False
    for iter in range(0, max_iter):
        drag_coef = calc_drag_coef(vel_guess, obj_dia, flu_kine_vis)
        vel_next = calc_obj_vel(obj_dia, obj_dens, flu_dens, drag_coef)
        if verbose:
            print 'Iter %d, Velocity %f m/s' % (iter+1, vel_next)
        if abs(vel_next-vel_guess) <= vel_thresh:
            vel_term = vel_next
            found = True
            break
        vel_guess = vel_next

    if not found:
        raise Exception('Failed to calculate the terminal velocity')

    if verbose:
        drag_force = calc_drag_force(drag_coef, obj_area, flu_dens, vel_term)
        force_bal = abs(obj_weight-buoy_force-drag_force)
        print 'Force balance: %g N' % force_bal
        print 'Terminal velocity: %f m/s' % vel_term

    return vel_term

def calc_obj_term_vel_2(time_res, force_bal_thresh, obj_dia, obj_dens,
        flu_dens, flu_kine_vis, verbose):
    # not a function of height
    obj_vol = pi*obj_dia**3/6
    obj_area = pi*obj_dia**2/4
    obj_mass = obj_dens*obj_vol
    obj_weight = obj_mass*grav
    buoy_force = flu_dens*grav*obj_vol

    time = 0
    fall = 0
    vel_term = 0
    force_bal = force_bal_thresh+1
    while force_bal > force_bal_thresh:
        time += time_res
        vel_grav = vel_term+grav*time_res
        drag_coef = calc_drag_coef(vel_grav, obj_dia, flu_kine_vis)
        drag_force = calc_drag_force(drag_coef, obj_area, flu_dens, vel_grav)
        force_bal = obj_weight-buoy_force-drag_force
        acc = force_bal/obj_mass
        vel_term += acc*time_res
        fall += vel_term*time_res

    if verbose:
        print 'Time: %f s' % time
        print 'Fall: %f m' % fall
        print 'Force balance: %g N' % force_bal
        print 'Terminal velocity: %f m/s' % vel_term

    return vel_term

print '=== Method 1 ==='
vel_term_1 = calc_obj_term_vel_1(vel_guess, vel_thresh, max_iter, obj_dia,
        obj_dens, flu_dens, flu_kine_vis, True)
print
print '=== Method 2 ==='
vel_term_2 = calc_obj_term_vel_2(time_res, force_bal_thresh, obj_dia, obj_dens,
        flu_dens, flu_kine_vis, True)
