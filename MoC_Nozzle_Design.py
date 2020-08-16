# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 12:00:30 2020

@author: ekrem
"""


from numericmethods import *
from math import *
import matplotlib.pyplot as plt


k=1.4
Me = 2.4
T0 = 486 #K
T_star = 404.2 #K
R = 287
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

def temp(Ma):
    return T0/(1+(k-1)/2*Ma**2)

##### DICTIONARY INITIALIZATION

point={}

theta_max = PM(Me)/2
Rt=0.1014640057401272                    #Throat Radius
y_axis=0


wall_points = [[y_axis, Rt]]

n = 24 # number of characteristic lines
delta_theta = theta_max -int(theta_max)
dt = int(theta_max)/(n)


theta_max_array = []

for h in range(n+1):
    theta_max_array.append(theta_max-h*dt)
theta_max_array.append(delta_theta)


changer=1
ind=0
indicator=n+1
old_wall_x = y_axis
old_wall_y = Rt

plt.figure(figsize=(10,6))

#%%%%%%%%%%%%%%%%%% SETTINGS  %%%%%%%%%%%%%%%%%

animation = False
figure_saving = False
print_wall = False
contour_text = False
temp_text = True

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    point[ind]={"Theta":theta_1, "Nu":nu_1, "Ma":Ma_1, "Mu":mu_1, "K-":K_1_minus, "K+":K_1_plus, "x":x1, "y":y1, "T":temp(Ma_1)}
    
    if animation==True:
        plt.pause(0.000001)

    plt.plot([x1, y_axis], [y_axis, Rt],'orange')


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
            
        plt.plot( [new_x, point[ind]['x']], [new_y,point[ind]['y']],'crimson')
        if animation==True:
            plt.pause(0.000001)
        ind+=1
        point[ind]={"Theta":new_theta, "Nu":new_nu,"Ma":new_Ma, "Mu":new_mu, "K-":right, "K+":left, "x":new_x, "y":new_y, "T":temp(new_Ma)}
  
    
    # WALL POINT
    if i == 0:
        alfa_wall_right = 0.5*(theta_max+point[ind]['Theta'])
        m_wall_right = tan(radians(alfa_wall_right))
    else:
        alfa_wall_right = 0.5*(point[mid_indicator]['Theta']+point[ind]['Theta'])
        m_wall_right = tan(radians(alfa_wall_right))
    
    alfa_wall_left = (point[ind]['Mu']+point[ind]['Theta'])
    m_wall_left = tan(radians(alfa_wall_left))

    wall_x = ((wall_points[-1][1]-point[ind]['y'])+m_wall_left*point[ind]['x']-m_wall_right*wall_points[-1][0])/(m_wall_left-m_wall_right)
    wall_y = point[ind]['y']+m_wall_left*(wall_x-point[ind]['x']) 

    if animation==True:
        plt.pause(0.000001)

    plt.plot( [wall_x, point[ind]['x']], [wall_y,point[ind]['y']],'forestgreen')

    
    ind+=1
    point[ind]={"Theta":point[ind-1]['Theta'], "Nu":point[ind-1]['Nu'],"Ma":point[ind-1]['Ma'], "Mu":point[ind-1]['Mu'], "K-":right, "K+":left, "x":wall_x, "y":wall_y, "T":temp(point[ind-1]['Ma'])}

    if animation==True:
        plt.pause(0.000001)
    
    if i==n:
        plt.plot( [wall_x, wall_points[-1][0]], [wall_y, wall_points[-1][1]],'indigo', linewidth=4, label='Wall Contour')
    else:
        plt.plot( [wall_x, wall_points[-1][0]], [wall_y, wall_points[-1][1]],'indigo', linewidth=4)

    
    wall_points.append([wall_x, wall_y])
       
    changer+=1
    indicator-=1
    
    print("Wall points: ",wall_x, "\t", wall_y,"\n")


plt.xlabel("Nozzle length (m) ",fontsize=16) #$\frac{x}{x0}$
plt.ylabel("Nozzle radius (m)",fontsize=16) #$\frac{y}{y0}$
# plt.title("Minimum Length Supersonic Nozzle",fontsize=20)

x_points = np.array([i[0] for i in wall_points])
y_points = np.array([i[1] for i in wall_points])

#''


###  REFLECTION

# centerline,=plt.plot([y_axis, x_points[-1]],[y_axis, y_axis],linewidth=1,color='gray')
# dashes=[50,5,3,5]
# centerline.set_dashes(dashes)
#plt.plot((x_points),(-y_points),'m')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plt.legend(loc=7, prop={'size': 12})

plt.show()

if figure_saving == True:
    file_name = "char"+str(n+1)+".pdf"
    plt.savefig(("plots/"+file_name), format='pdf',dpi=1200)    
    print("Image file ",file_name," has been generated!")
   
if print_wall == True:
    
    zw = [0]*len(wall_points)

    # UNIT CONVERSION
    
    #for m to mm -> mm=1000
    #    m to m  -> mm=1
    
    mm=1 ##m to mm
    
    with open("pointss.asc", "w") as o:
        count  = 0
        for i in range(len(wall_points)):
            
            if i ==len(wall_points):
                print(wall_points[i][0]*mm," ", wall_points[i][1]*mm," ", zw[i]*mm, file=o)
            else:
                print(wall_points[i][0]*mm," ", wall_points[i][1]*mm,  zw[i]*mm, file=o,end = "\n")
            count+=1

if contour_text == True:
    file_name = "machcontour"+str(n +1) + ".csv"
    with open(("point_folder/"+file_name), "w") as m:
        z_axis = 0
        print('x' , 'y', 'Ma', file=m)
        print(y_axis, Rt, 1, file=m)
        
        
        
        for i in point:

            if i == len(point):
                print(point[i]['x'] , point[i]['y'], point[i]['Ma'], file=m)
            else:
                print(point[i]['x'] , point[i]['y'], point[i]['Ma'], file=m, end="\n")
        
        def last_line(x):
            return point[len(point)]['y']/(point[len(point)]['x']-point[len(point)-1]['x'])*(-point[len(point)-1]['x']+x)
        arr = np.linspace(point[len(point)-1]['x'], point[len(point)]['x'],500)
        for i in arr:    
            print(i, last_line(i), Me,  file=m)
        print(x_points[-1], 0, Me, file=m)
        # print(0, 0, 0.9, file=m)
        print("---> Mach contour text file saved as ",file_name)

if temp_text == True:
    file_name = "tempcontour"+str(n +1) + ".csv"
    with open(("point_folder/"+file_name), "w") as m:
        z_axis = 0
        print('x' , 'y', 'T', file=m)
        print(y_axis, Rt, T_star, file=m)
        
        
        
        for i in point:

            if i == len(point):
                print(point[i]['x'] , point[i]['y'], point[i]['T'], file=m)
            else:
                print(point[i]['x'] , point[i]['y'], point[i]['T'], file=m, end="\n")
        
        def last_line(x):
            return point[len(point)]['y']/(point[len(point)]['x']-point[len(point)-1]['x'])*(-point[len(point)-1]['x']+x)
        arr = np.linspace(point[len(point)-1]['x'], point[len(point)]['x'],500)
        for i in arr:    
            print(i, last_line(i), 218.4014869888476,  file=m)
        print(x_points[-1], 0, 218.4014869888476, file=m)
        # print(0, 0, 0.9, file=m)
        print("---> Temperature contour text file saved as ",file_name)
