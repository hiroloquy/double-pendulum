# Setting ----------------------------------------
reset
set nokey
set term gif animate delay 4 size 854,480
set output 'dp_trajectory.gif'

set xr[-40:40]
set yr[-40:40]
set xl 'x' font 'Times New Roman:Italic, 20'
set yl 'y' font 'Times New Roman:Italic, 20'
set tics font 'Times New Roman,18'

set size ratio -1
set grid

# Parameter ----------------------------------------
m1  = 1.0              # mass of bob
m2  = 1.5
l1  = 23.0             # length of link
l2  = 16.0
g  = 9.81              # gravitational acceleration
dt   = 0.01            # step
M  = m2 / (m1+m2)      # coefficient
l  = l2 / l1
omega2 = g / l1
dx  = dt/6
limit = 4000           # loop limit
R  = 1                 # radius of bob
r  = 0.1               # radius of trajectory

# Runge-Kutta 4th ----------------------------------------
f1(a, b, c, d) = c     # theta1'
f2(a, b, c, d) = d     # theta2'
f3(a, b, c, d) = \     # theta1''
 (omega2*l*(-sin(a)+M*cos(a-b)*sin(b))-M*l*(c**2*cos(a-b)+l*d**2)*sin(a-b)) / (l-M*l*(cos(a-b))**2)
f4(a, b, c, d) = \     # theta2''
 (omega2*cos(a-b)*sin(a)-omega2*sin(b)+(c**2+M*l*d**2*cos(a-b))*sin(a-b)) / (l-M*l*(cos(a-b))**2)

# Function using label ----------------------------------------
#Parameter
label(a, b, c, d, e, f, g, h) = sprintf("\
{/Times:Italic m_{/Times:Normal 1}} = %.2f [kg]\n\
{/Times:Italic m_{/Times:Normal 2}} = %.2f [kg]\n\
{/Times:Italic l_{/Times:Normal 1}}   = %3.1f [m]\n\
{/Times:Italic l_{/Times:Normal 2}}   = %3.1f [m]\n\
{/Times:Italic g}    = %3.2f [m/s^2]\n\
{/Times:Italic dt}   = %.2f [s]\n\n\
{/symbol-oblique q_{/Times:Normal 01}}  = %3.2f [rad]\n{/symbol-oblique q_{/Times:Normal 02}}  = %3.2f [rad]", a, b, c, d, e, f, g, h) 

# Time
time(t) = sprintf("{/Times:Italic t} = %5.2f [s]", t)

# Plot ----------------------------------------
# Initial Value
x1 = 3.5*pi/6         # theta1
x2 = 4.5*pi/6         # theta2
x3 = 0.0              # theta1'
x4 = 0.0              # theta2'
t  = 0.0              # time

# Draw initiate state for 70 steps
do for [i = 1:70] {
    # Links
    set arrow 1 nohead lw 2 from 0, 0 to l1*sin(x1), -l1*cos(x1) lc -1
    set arrow 2 nohead lw 2 from l1*sin(x1), -l1*cos(x1) to l1*sin(x1)+l2*sin(x2), -(l1*cos(x1)+l2*cos(x2)) lc -1

    # Time and Parameter
    set title time(t) font 'Times:Normal, 20'
    set label 1 left at 45, 20 label(m1, m2, l1, l2, g, dt, x1, x2) font 'Times:Normal, 18'

    # Anchor
    set object 1 circle at 0, 0 fc rgb 'black' size R fs solid front

    # Bobs
    set object 2 circle at l1*sin(x1), -l1*cos(x1) fc rgb 'blue' size R fs solid  front
    set object 3 circle at l1*sin(x1)+l2*sin(x2), -(l1*cos(x1)+l2*cos(x2)) fc rgb 'red' size R fs solid front

    # Draw
    plot 1/0

    # Trajectory turns smaller
    set object 2 size r
    set object 3 size r
}

# Draw simulation for limit steps
do for [i = 1:limit] {
    # Calculate using Runge-Kutta 4th
    t = t + dt
    k11 = f1(x1, x2, x3, x4)
    k12 = f2(x1, x2, x3, x4)
    k13 = f3(x1, x2, x3, x4)
    k14 = f4(x1, x2, x3, x4)
    k21 = f1(x1+dt/2*k11, x2+dt/2*k12, x3+dt/2*k13, x4+dt/2*k14 )
    k22 = f2(x1+dt/2*k11, x2+dt/2*k12, x3+dt/2*k13, x4+dt/2*k14 )
    k23 = f3(x1+dt/2*k11, x2+dt/2*k12, x3+dt/2*k13, x4+dt/2*k14 )
    k24 = f4(x1+dt/2*k11, x2+dt/2*k12, x3+dt/2*k13, x4+dt/2*k14 )
    k31 = f1(x1+dt/2*k21, x2+dt/2*k22, x3+dt/2*k23, x4+dt/2*k24 )
    k32 = f2(x1+dt/2*k21, x2+dt/2*k22, x3+dt/2*k23, x4+dt/2*k24 )
    k33 = f3(x1+dt/2*k21, x2+dt/2*k22, x3+dt/2*k23, x4+dt/2*k24 )
    k34 = f4(x1+dt/2*k21, x2+dt/2*k22, x3+dt/2*k23, x4+dt/2*k24 )
    k41 = f1(x1+dt*k31,  x2+dt*k32,  x3+dt*k33,  x4+dt*k34 )
    k42 = f2(x1+dt*k31,  x2+dt*k32,  x3+dt*k33,  x4+dt*k34 )
    k43 = f3(x1+dt*k31,  x2+dt*k32,  x3+dt*k33,  x4+dt*k34 )
    k44 = f4(x1+dt*k31,  x2+dt*k32,  x3+dt*k33,  x4+dt*k34 )

    # Update angle and angular velocity
    x1 = x1 + dx * (k11 + 2*k21 + 2*k31 + k41)
    x2 = x2 + dx * (k12 + 2*k22 + 2*k32 + k42)
    x3 = x3 + dx * (k13 + 2*k23 + 2*k33 + k43)
    x4 = x4 + dx * (k14 + 2*k24 + 2*k34 + k44)

    # Draw links
    set arrow 1 nohead lw 2 from 0, 0 to l1*sin(x1), -l1*cos(x1) lc -1
    set arrow 2 nohead lw 2 from l1*sin(x1), -l1*cos(x1) to l1*sin(x1)+l2*sin(x2), -(l1*cos(x1)+l2*cos(x2)) lc -1

    # Draw bobs
    set object 2*i+2 circle at l1*sin(x1), -l1*cos(x1) fc rgb "blue" size R fs solid front
    set object 2*i+3 circle at l1*sin(x1)+l2*sin(x2), -(l1*cos(x1)+l2*cos(x2)) fc rgb "red" size R fs solid front

    # Draw time
    set title time(t)

    # Decimate and draw
    if (i%20==0){
        plot 1/0
    }

    # Trajectory turns smaller
    set object 2*i+2 size r
    set object 2*i+3 size r
}

set out