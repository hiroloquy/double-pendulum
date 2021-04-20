# Setting ----------------------------------------
reset
set term gif animate delay 4 size 854,480
set output 'dp_afterimage.gif'

set nokey
set xr[-40:40]
set yr[-40:40]
set xl 'x' font 'Times New Roman:Italic, 20'
set yl 'y' font 'Times New Roman:Italic, 20'
set tics font 'Times New Roman,18'

set size ratio -1
set grid

# Parameter ----------------------------------------
m1  = 1.0               # mass
m2  = 1.5
l1  = 23.0              # length of link
l2  = 16.0
g  = 9.81               # gravitational acceleration
M  = m2 / (m1+m2)       # coefficient
l  = l2 / l1
omega2 = g / l1
dt   = 0.08             # step
limit = 2000            # limit of loop
cnt  = 10               # nunber of trajectory
dis  = 5                # start to disappear
r  = 1                  # radius

# Runge-Kutta 4th ----------------------------------------
f1(a, b, c, d) = c      # theta1'
f2(a, b, c, d) = d      # theta2'
f3(a, b, c, d) = \      # theta1''
 (omega2*l*(-sin(a)+M*cos(a-b)*sin(b))-M*l*(c**2*cos(a-b)+l*d**2)*sin(a-b)) / (l-M*l*(cos(a-b))**2)
f4(a, b, c, d) = \      # theta2''
 (omega2*cos(a-b)*sin(a)-omega2*sin(b)+(c**2+M*l*d**2*cos(a-b))*sin(a-b)) / (l-M*l*(cos(a-b))**2)

# Function of Text ----------------------------------------
#Parameter
label(a, b, c, d, e, f, g, h) = sprintf("\
{/Times:Italic m_{/Times:Normal 1}} = %.2f [kg]\n\
{/Times:Italic m_{/Times:Normal 2}} = %.2f [kg]\n\
{/Times:Italic l_{/Times:Normal 1}}   = %3.1f [m]\n\
{/Times:Italic l_{/Times:Normal 2}}   = %3.1f [m]\n\
{/Times:Italic g}    = %3.2f [m/s^2]\n\
{/Times:Italic dt}   = %.2f [s]\n\n\
{/symbol-oblique q_{/Times:Normal 01}}  = %3.2f [rad]\n\
{/symbol-oblique q_{/Times:Normal 02}}  = %3.2f [rad]", a, b, c, d, e, f, g, h) 

#Time
time(t) = sprintf("{/Times:Italic t} = %5.2f [s]", t)

# Plot ----------------------------------------
# Initial Value
x1 = 3*pi/6             # theta1
x2 = 4*pi/6             # theta2
x3 = 0.0                # theta1'
x4 = 0.0                # theta2'
t  = 0.0                # time

# Draw initiate state for 70 steps
do for [i = 1:70] {
    # Links
    set arrow 1 nohead lw 2 from 0, 0 to l1*sin(x1), -l1*cos(x1) lc -1
    set arrow 2 nohead lw 2 from l1*sin(x1), -l1*cos(x1) to l1*sin(x1)+l2*sin(x2), -(l1*cos(x1)+l2*cos(x2)) lc -1

    # Time and Parameter
    set title time(t) font 'Times:Normal, 20'
    set label 1 left at 45, 20 label(m1, m2, l1, l2, g, dt, x1, x2) font 'Times:Normal, 18'

    # Anchor
    set object 1 circle at 0, 0 fc rgb 'black' size r fs solid front

    # Bobs
    set object 2 circle at l1*sin(x1), -l1*cos(x1) fc rgb 'blue' size r fs solid  front
    set object 3 circle at l1*sin(x1)+l2*sin(x2), -(l1*cos(x1)+l2*cos(x2)) fc rgb 'red' size r fs solid front

    # Draw
    plot 1/0
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

    # Update theta
    x1 = x1 + dt/6*(k11+2*k21+2*k31+k41)
    x2 = x2 + dt/6*(k12+2*k22+2*k32+k42)
    x3 = x3 + dt/6*(k13+2*k23+2*k33+k43)
    x4 = x4 + dt/6*(k14+2*k24+2*k34+k44)

    # Update links
    set arrow 1 nohead lw 2 from 0, 0 to l1*sin(x1), -l1*cos(x1) lc -1
    set arrow 2 nohead lw 2 from l1*sin(x1), -l1*cos(x1) to l1*sin(x1)+l2*sin(x2), -(l1*cos(x1)+l2*cos(x2)) lc -1

    # Draw bobs
    set object 2*i+2 circle at l1*sin(x1), -l1*cos(x1) fc rgb "blue" size r fs solid front
    set object 2*i+3 circle at l1*sin(x1)+l2*sin(x2), -(l1*cos(x1)+l2*cos(x2)) fc rgb "red" size r fs solid front

    # Update time
    set title time(t)

    # Start to disappear
    if(i>=dis){
        lim2 = cnt - dis - 1             # Number of trajectory turning lighter

        do for [j=1:lim2] {
            k1 = 2*(i-dis+j)
            k2 = k1+1

            tp = 1.0 - 0.15*(lim2-j+1)   # Decrement transparent and size of trajectory by 0.15
            si = r - 0.15*(lim2-j+1)

            set object k1 size si fs transparent solid tp noborder
            set object k2 size si fs transparent solid tp noborder
        }
    }

    if(i>=cnt){                       # Remove old objects
     k1 = 2*(i-cnt)+2
     k2 = k1+1
     unset object k1
     unset object k2
    }

    plot 1/0                             # Draw
}

set out