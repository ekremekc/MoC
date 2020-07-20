# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 12:00:30 2020

@author: ekrem
"""


from numericmethods import *
from math import *
import matplotlib.pyplot as plt
import matplotlib

k=1.4
Me = 2.4

spec_heat = k

def PM(M):
    ratio = (spec_heat+1)/(spec_heat-1)
    return sqrt(ratio)*degrees(atan(sqrt(1/ratio*(M**2-1))))-degrees(atan(sqrt(M**2-1)))

def ma_finder(pm_angle):
    def zero(M):
        return PM(M)-pm_angle
    return bisection_method(zero, 1, 10)

def mu(M):
    return degrees(asin(1/M))


##### ARRAYS INITIALIZATION

point={}



theta_max = PM(Me)/2
Rt=1
y_axis=0

wall_points = [[y_axis, Rt]]

n = 100 # number of characteristic lines
delta_theta = theta_max -int(theta_max)
dt = int(theta_max)/(n)


theta_max_array = []

for h in range(n+1):
    theta_max_array.append(theta_max-h*dt)
theta_max_array.append(delta_theta)
# theta_max_array = np.linspace(theta_max, delta_theta,n+1)

changer=1
ind=0
indicator=n+1
old_wall_x = y_axis
old_wall_y = Rt

plt.figure(figsize=(16,12))

animation = False

if animation==True:
    plt.plot(1)
    plt.show(block=False)
    plt.pause(4)

for i in range(0,n+1):
    
    # AXIS POINT
    theta_a = delta_theta+i*dt
    nu_a = PM(1)+theta_a
    Ma_a = ma_finder(nu_a)
    mu_a = mu(Ma_a)

    K_a_minus = theta_a + nu_a
    K_1_minus = K_a_minus
    K_1_plus = -theta_a - nu_a 
    
    theta_1 = 0.5*(K_1_minus+K_1_plus)
    
    nu_1 = 0.5*(K_1_minus-K_1_plus)
    Ma_1 = ma_finder(nu_1)
    mu_1 = mu(Ma_1)
      
    if i ==0:
        alfa_a_1 = 0.5*(theta_a+theta_1)-0.5*(mu_a+mu_1)
    else:
        alfa_a_1 = 0.5*(point[ind-indicator]["Theta"]+theta_1)-0.5*(point[ind-indicator]["Mu"]+mu_1)
    
    m_a_1 = tan(radians(alfa_a_1))
    
    if i==0:
        x1 = -Rt/m_a_1
        y1 = y_axis
        print("Axis Points: ",x1,y1)
    else:
        x1=point[ind-indicator]['x']-point[ind-indicator]['y']/m_a_1
        y1 = y_axis
        print("Axis Points: ",x1,y1)
    
    ind+=1
    point[ind]={"Theta":theta_1, "Nu":nu_1, "Ma":Ma_1, "Mu":mu_1, "K-":K_1_minus, "K+":K_1_plus, "x":x1, "y":y1}
    
    if animation==True:
        plt.pause(0.000001)

    plt.plot([x1, y_axis], [y_axis, Rt],'coral')


    mid_indicator = ind-indicator
    
    #INTERNAL POINTS
    for j in range(changer,n+1):
        theta = delta_theta+j*dt
        nu = PM(1)+theta
        Ma = ma_finder(nu)
        mu_ = mu(Ma)
        right = theta + nu
        left = K_1_plus
        
        new_theta = 0.5*(right + left)
        new_nu    = 0.5*(right - left) 
        new_Ma    = ma_finder(new_nu)
        new_mu    = mu(new_Ma)
         
        new_alfa_right = 0.5*(theta+new_theta)-0.5*(mu_+new_mu)
        m_alfa_right = tan(radians(new_alfa_right))
        
        new_alfa_left  = 0.5*(point[ind]['Theta']+new_theta)+0.5*(point[ind]['Mu']+new_mu) 
        m_alfa_left = tan(radians(new_alfa_left))
        
        if i==0:
            new_x = ((Rt-point[ind]['y'])+m_alfa_left*point[ind]['x']-m_alfa_right*y_axis)/(m_alfa_left-m_alfa_right)
            new_y = point[ind]['y']+m_alfa_left*(new_x-point[ind]['x'])
        else:
            new_x = ((point[mid_indicator]['y']-point[ind]['y'])+m_alfa_left*point[ind]['x']-m_alfa_right*point[mid_indicator]['x'])/(m_alfa_left-m_alfa_right)
            new_y = point[ind]['y']+m_alfa_left*(new_x-point[ind]['x'])
            mid_indicator+=1
            
        plt.plot( [new_x, point[ind]['x']], [new_y,point[ind]['y']],'deepskyblue')
        if animation==True:
            plt.pause(0.000001)
        ind+=1
        point[ind]={"Theta":new_theta, "Nu":new_nu,"Ma":new_Ma, "Mu":new_mu, "K-":right, "K+":left, "x":new_x, "y":new_y}

        # print(new_x,"\t", new_y)    
    
    # WALL POINT
    if i == 0:
        alfa_wall_right = 0.5*(theta_max+point[ind]['Theta'])
        m_wall_right = tan(radians(alfa_wall_right))
    else:
        # print("MID INDICATOR: ",mid_indicator, point[mid_indicator]['Theta'] )
        alfa_wall_right = 0.5*(point[mid_indicator]['Theta']+point[ind]['Theta'])
        m_wall_right = tan(radians(alfa_wall_right))
    
    alfa_wall_left = (point[ind]['Mu']+point[ind]['Theta'])
    m_wall_left = tan(radians(alfa_wall_left))

    wall_x = ((wall_points[-1][1]-point[ind]['y'])+m_wall_left*point[ind]['x']-m_wall_right*wall_points[-1][0])/(m_wall_left-m_wall_right)
    wall_y = point[ind]['y']+m_wall_left*(wall_x-point[ind]['x']) 

    if animation==True:
        plt.pause(0.000001)

    plt.plot( [wall_x, point[ind]['x']], [wall_y,point[ind]['y']],'olive')

    
    ind+=1
    point[ind]={"Theta":point[ind-1]['Theta'], "Nu":point[ind-1]['Nu'],"Ma":point[ind-1]['Ma'], "Mu":point[ind-1]['Mu'], "K-":right, "K+":left, "x":wall_x, "y":wall_y}

    if animation==True:
        plt.pause(0.000001)
    
    plt.plot( [wall_x, wall_points[-1][0]], [wall_y, wall_points[-1][1]],'m')

    
    wall_points.append([wall_x, wall_y])
       
    changer+=1
    indicator-=1
    
    print("Wall points: ",wall_x, "\t", wall_y,"\n")

# plt.ylim(0,wall_points[-1][0]/2)
plt.xlabel("Nozzle Length")
plt.ylabel("Nozzle Radius")
plt.title("Minimum Length Nozzle")

x_points = np.array([i[0] for i in wall_points])
y_points = np.array([i[1] for i in wall_points])

centerline,=plt.plot([y_axis, x_points[-1]],[y_axis, y_axis],linewidth=1,color='gray')
dashes=[50,5,3,5]
centerline.set_dashes(dashes)


plt.plot((x_points),(-y_points),'m')

plt.show()


# zw = [0]*len(wall_points)

# # UNIT CONVERSION

# #for m to mm -> mm=1000
# #    m to m  -> mm=1

# mm=1 ##m to mm

# with open("pointss.asc", "w") as o:
#     count  = 0
#     for i in range(len(wall_points)):
        
#         if i ==len(wall_points):
#             print(wall_points[i][0]*mm," ", wall_points[i][1]*mm," ", zw[i]*mm, file=o)
#         else:
#             print(wall_points[i][0]*mm," ", wall_points[i][1]*mm,  zw[i]*mm, file=o,end = "\n")
#         count+=1


