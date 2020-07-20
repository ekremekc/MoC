import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import math


""" 

      #############          GENERAL GAS CLASS         ###########

"""

class Gas():
    def __init__(self,T,P,R_u=8.31447):

        self.T = T
        self.P = P
        self.R_u=R_u
        self.normalshock=self.Shock(self)
        
    def gas_list(self):
        print(" Code\t","Gas","\n",
              "----\t","---","\n",
              "Air\t","Air" ,"\n",
              "Ar\t\t","Argon" ,"\n" ,
              "CO2\t","Carbon dioxide" ,"\n",
              "CO\t\t","Carbon monoxide" ,"\n",
              "N2\t\t","Nitrogen" ,"\n"   
                
              )
        
    def area(self, diameter):
        return (np.pi*(diameter**2))/4
    def critical_area(self,massflowrate):
        return massflowrate/(self.P*1000*(self.k**(1/2))*(2/(self.k+1))**((self.k+1)/(2*self.k-2))/((self.R*1000*self.T)**(1/2)))
    def critical_m_dot(self, Ma, diameter=1):
        return self.critical_density()*self.area(diameter)*self.critical_speed_of_sound(Ma)
    def critical_temperature(self, Ma):
        return self.stagnation_temp(Ma)*2/(self.k+1)
    def critical_pressure(self):
        return self.P*(2/(self.k+1))**(self.k/(self.k-1))
    def critical_density(self):
        return self.rho*(2/(self.k+1))**(1/(self.k-1))
    def critical_speed_of_sound(self, Ma):
        return np.sqrt(self.k*self.R*self.critical_temperature(Ma)*1000)
    def density(self):
        return self.P/(self.R*self.T)
    def diameter(self, area):
        return np.sqrt(4/np.pi*area)
    def enthalpy(self):
        return self.cp*self.T
    def exit_temperature(self,Mach):
        return self.T/(1+(self.k-1)/2*Mach**2)
    def exit_pressure(self,Mach):
        return self.P/(1+(self.k-1)/2*Mach**2)**(self.k/(self.k-1))
    def exit_density(self, Mach):
        return self.rho/(1+(self.k-1)/2*Mach**2)**(1/(self.k-1))
    def exit_speed(self, Mach):
        return Mach*np.sqrt(self.k*self.R*self.exit_temperature(Mach)*1000)
    def exit_area(self, Throat_Area, Mach):
        return Throat_Area*(1/Mach)*((2/(self.k+1))*(1+(self.k-1)/2*Mach**2))**((self.k+1)/(2*self.k-2))
    def mach_number(self, velocity):
        return velocity/self.speed_of_sound()
    def m_dot(self, velocity, diameter=1):
        return self.density()*self.area(diameter)*velocity   
    def mfr(self,velocity, diameter):
        return self.critical_pressure()*self.area(diameter)*self.mach_number(velocity)*np.sqrt(self.k/(self.R*self.critical_temperature())) 
    def mass_flowrate(self, velocity, diameter=1):
        return (self.area(diameter)*self.mach_number(velocity)*self.stagnation_pressure(velocity)*np.sqrt(self.k/(self.R*self.stagnation_temp(velocity))))\
              /((1+(self.k-1)*(self.mach_number(velocity)**2)/2)**((self.k+1)/(2*(self.k-1))))

    def ma_finder(self, section,  area_ratio, show_iterations=False, tolerance=10e-6, method="bisection"):
        try:
            if section !="upward" and section !="downward":
                raise NameError("Please specify the flow by using these keywords:  \"upward\" or \"downward\"")
            def finder(Ma):
                value = (1/Ma*((1+0.5*(self.k-1)*Ma**2)/(0.5*(self.k+1)))**(0.5*(self.k+1)/(self.k-1)))
                if method=='golden' or method=='secant':
                    target = abs(value - area_ratio)
                elif method=='bisection':
                    target = value - area_ratio
                return target
            # def check_boundaries(Ma_0, Ma_1):
            #     if section=="upward":
            #         if Ma_0>1 or Ma_1>1:
            #             Ma_0 = 1/Ma_0
            #             Ma_1 = Ma_0+0.001
            #             # print("ma kucuk 1 den calisti")
            #     elif section=="downward":
            #         if Ma_0<1 or Ma_1<1:
            #             Ma_0 = 1+Ma_0
            #             Ma_1 = Ma_0+0.1
            #             # print("ma buyuk 1 den calisti")
                
            
            if section=="upward":
                if method=='bisection':
                    Ma=bisection_method( finder,0, 1, tolerance = 10e-6,show_iterations=show_iterations)
                elif method=='secant':
                    Ma=secant_method( finder,0, 1, tolerance = 10e-6,show_iterations=show_iterations)
                elif method=='golden':
                    Ma=golden_section(finder,0, 1, tolerance = 10e-6,show_iterations=show_iterations)
                
             
            elif section=="downward":
                if method=='bisection':
                    Ma=bisection_method( finder,1, 5, tolerance = 10e-6,show_iterations=show_iterations)
                elif method=='secant':
                    Ma=secant_method( finder,1, 5, tolerance = 10e-6,show_iterations=show_iterations)
                elif method=='golden':
                    Ma=golden_section(finder,1, 5, tolerance = 10e-6,show_iterations=show_iterations)
            return Ma
    
        except NameError:       
            raise NameError("Please specify the flow by using these keywords:  \"upward\" or \"downward\"") from None
        except ValueError:       
            raise ValueError("Given area is smaller than throat area. Program has terminated.\n Hint: You could change the division number.") from None
    
    def plot(self,area_start, area_end, Mach_start, y_axis='T', color_bar='Ma', division=250 ,x_axis='A', method="bisection"):
        area_upward = np.linspace(area_start, self.throat_area(area_start,Mach_start), division)
        area_downward = np.linspace(self.throat_area(area_start,Mach_start), area_end, division)
        area_total = np.concatenate((area_upward,area_downward))

        
        ST = self.stagnation_temp(Mach_start)
        temp_upward = []
        Ma_upward = []
        
        for i in range(division):
            ratio = self.throat_area_ratio(area_upward[i], area_start, Mach_start)
            
            Ma=self.ma_finder("upward",ratio,method=method)

            Ma_upward.append(Ma)
                      
            temp_upward.append(self.temperature(Ma, ST))
        
        
        temp_downward = []
        Ma_downward = []
        
        for i in range(division):
            ratio = self.throat_area_ratio(area_downward[i], area_start, Mach_start)
           
            Ma=self.ma_finder("downward",ratio,method=method)
            Ma_downward.append(Ma)
            temp_downward.append(self.temperature(Ma, ST))

        temp_total = temp_upward +temp_downward
        Ma_total = Ma_upward +Ma_downward
        
        fig = plt.figure(figsize=(10,7.5))
        ax = fig.add_subplot(111)
        xs = np.linspace(0,1,2*division)

        
        if y_axis == 'T':
            y_lbl='Temperature (K)'
            if color_bar=='Ma':
                color = Ma_total
                mp = ax.scatter((xs),(temp_total),c=color,cmap=plt.cm.get_cmap('jet'))
                c_lbl = 'Mach Number'
            elif color_bar=='T':
                mp = ax.scatter((xs),(temp_total),c=temp_total,cmap=plt.cm.get_cmap('jet'))
                c_lbl = 'T (K)'
        elif y_axis == 'Ma':
            y_lbl='Mach Number'
            if color_bar=='Ma':
                color = Ma_total
                mp = ax.scatter((xs),(Ma_total),c=color,cmap=plt.cm.get_cmap('jet'))
                c_lbl = 'Mach Number'
            elif color_bar=='T':
                mp = ax.scatter((xs),(Ma_total),c=temp_total,cmap=plt.cm.get_cmap('jet'))
                c_lbl = 'T (K)'
                
        cb = plt.colorbar(mp)
        cb.set_label(c_lbl)
        ax.set(title=r'Converging- Diverging Nozzle',
                       xlabel='Area $m^2$', ylabel=y_lbl)
            
    
        tick_labels=[]
        
        for j in np.linspace(0,(2*division),7):
            if j==2*division:
                tick_labels.append(round(area_total[-1],4))
            else:    
                tick_labels.append(round(area_total[int(j)],4))
        plt.xticks(np.linspace(0,1,7),tick_labels)
        
        
        plt.show()


    def pressure(self, Mach, Stagnation_Pressure):
        return Stagnation_Pressure/((1+0.5*(self.k-1)*Mach**2)**(self.k/(self.k-1)))
    def speed_of_sound(self):
        return np.sqrt(self.k*self.R*self.T*1000)
    def stagnation_temp(self,Mach):
        return self.T*(1+(self.k-1)/2*Mach**2)
    def stagnation_pressure(self,Mach):
        return self.P*(1+0.5*(self.k-1)*Mach**2)**(self.k/(self.k-1))
    def temperature(self, Mach, Stagnation_Temperature):
        return Stagnation_Temperature/(1+(self.k-1)/2*Mach**2)
    def throat_area(self,known_area,Mach):
        return known_area/((1/Mach)*((2/(self.k+1))*(1+(self.k-1)/2*Mach**2))**((self.k+1)/(2*self.k-2)))
    def throat_area_ratio(self,wanted_area, known_area,known_Mach):
        return wanted_area/self.throat_area(known_area, known_Mach)
    
    class Shock():
        def __init__(self, gas):
            self.gas = gas
            
                
        def P2(self, Ma1, P1):
            return P1*(1/(self.gas.k+1)*(2*self.gas.k*Ma1**2-(self.gas.k-1)))
        def Ma2(self,Ma1):
            return np.sqrt(((self.gas.k-1)*Ma1**2+2)/(2*self.gas.k*Ma1**2-(self.gas.k-1)))
        def P0_2(self,Stagnation_Pressure, Ma1):
            return Stagnation_Pressure*((((self.gas.k+1)*Ma1**2)/(2+(self.gas.k-1)*Ma1**2))**(self.gas.k/(self.gas.k-1))\
                                        *((self.gas.k+1)/(2*self.gas.k*Ma1**2-(self.gas.k-1)))**(1/(self.gas.k-1)))
        def area_shock_star(self, area1_star, Ma1):
            return area1_star*(self.Ma2(Ma1)/Ma1)*((2+(self.gas.k-1)*Ma1**2)/(2+(self.gas.k-1)*self.Ma2(Ma1)**2))**((self.gas.k+1)/(2*self.gas.k-2))
        
        def Ma_beforeshock(self, P2_P1):
            return np.sqrt((P2_P1*(self.gas.k+1)+(self.gas.k-1))/(2*self.gas.k))
        
        def T2(self,T1,Ma1):
            return T1*(2+(self.gas.k-1)*Ma1**2)*(2*self.gas.k*Ma1**2-(self.gas.k-1))/(((self.gas.k+1)**2)*(Ma1**2))
        def V2(self, T1, V1):
            return np.sqrt(2*self.gas.cp*(T1-self.T2(T1, V1/(self.gas.speed_of_sound())))+V1**2)
    

    
class Air(Gas):
    def __init__(self,T=298.15,P=101.325):
        super().__init__(T, P)  
        self.M=28.97
        self.k=1.4
        self.R=self.R_u/self.M
        self.cp=1.9327E-10*self.T**4 - 7.9999E-07*self.T**3 + 1.1407E-03*self.T**2 - 4.4890E-01*self.T + 1.0575E+03
        self.rho = self.P/(self.R*self.T) 

 
class CO2(Gas):
    def __init__(self,T=298.15,P=101.325):
        super().__init__(T, P)  
        self.M=44.01
        self.k=1.289
        self.R=self.R_u/self.M
        self.cp=0.849
        self.rho = self.P/(self.R*self.T)
        
class CO(Gas):
    def __init__(self,T=298.15,P=101.325):
        super().__init__(T, P)  
        self.M=28.01
        self.k=1.4
        self.R=self.R_u/self.M
        self.cp=1.039
        self.rho = self.P/(self.R*self.T)
 
class N2(Gas): 
    def __init__(self,T=298.15,P=101.325):
        super().__init__(T, P) 
        self.M=28.01
        self.k=1.4
        self.R=self.R_u/self.M
        self.cp=1.040
        self.rho = self.P/(self.R*self.T)

class Ar(Gas):
    def __init__(self,T=298.15,P=101.325):
        super().__init__(T, P) 
        self.M=39.95
        self.k=1.667
        self.R=self.R_u/self.M
        self.cp=0.5203
        self.rho = self.P/(self.R*self.T)

""" 

      #############          NUMERICAL METHODS         ###########

"""

def golden_section(func,starting, ending, show_iterations=False, tolerance = 10e-6):
    gr=(np.sqrt(5)+1)/2-1
    dm=tolerance
    a0 = starting+dm
    b0 = ending-dm
    count=0
    while True:
        count+=1
        # print(finder(Ma_0))
        # print(finder(Ma_1))
        d=gr*(b0-a0)
        a1=a0+d
        b1=b0-d
        if abs((a1-b1)/a1)<=tolerance:
            if 1>=ending:    
                print("The Mach number below unity is: ",a1,"\n")
            elif starting>=1:    
                print("The Mach number above unity is: ",a1,"\n")
            break
        else:
            if func(a0)>func(b0):
                a0=a1
                b1=b1
                
            else:
                a0=a0
                b0=b1

        if show_iterations ==True:
            print("Iteration ", count, " :",a1)
    return (a1+b1)/2

def secant_method(func, lower_bound, upper_bound, show_iterations=False,tolerance=10e-6):
    
    Ma_0 = (upper_bound+lower_bound)/2
    dMa = 0.01
    Ma_1 = Ma_0+dMa
    count=0
    while True:
        count+=1
        Ma_2 = Ma_1 - func(Ma_1)*(Ma_1-Ma_0)/(func(Ma_1)-func(Ma_0))
        if show_iterations ==True:
            print("Iteration ", count, " :",Ma_2)
        if func(Ma_2)<=tolerance:
            if show_iterations ==True:
                print("The Mach number below unity is: ",Ma_2,"\n")
            break
        else:
            Ma_0 = Ma_1
            Ma_1 = Ma_2 
    return Ma_2

def bisection_method(func, lower_bound, upper_bound, show_iterations=False,tolerance=10e-6):
    if lower_bound==0 :
        lower_bound+=tolerance
    a=lower_bound
    b= upper_bound
    count = 0
    while True:
        count+=1
        c = (a+b)/2
        if abs(func(c))<=tolerance:
            if show_iterations ==True:
                print("The root is: ",c,"\n")
            break
        else:
            if   func(a)*func(c)>func(b)*func(c):
                b=b
                a=c     
            else:     
                a=a  
                b=c
        if show_iterations ==True:
            print("Iteration ", count, " :",c)
    return c

""" 

      #############          ROCKET NOZZLE CLASS         ###########

"""
class Nozzle(Gas):
    def __init__(self, class_gas):
        self.T=class_gas.T
        self.P=class_gas.P
        self.k=class_gas.k
        self.M=class_gas.M
        self.k=class_gas.k
        self.R=class_gas.R_u/class_gas.M
        self.cp=class_gas.cp
        self.rho = class_gas.P/(class_gas.R*class_gas.T)
    def critical_throat_pressure(self):
        return self.P*(2/(self.k+1))**(self.k/(self.k-1))
    def exit_mach(self,backflow_pressure):
        if self.ischoked(backflow_pressure):
            Ma = 1
        else:
            Ma = np.sqrt(5*((self.P/backflow_pressure)**(2/7)-1))
        return Ma
    def ischoked(self, backflow_pressure ):
        if backflow_pressure < self.critical_pressure():
            condition=True
        else:
            condition = False
        return condition
    def massflowrate(self, backflow_pressure, area):
        if self.ischoked(backflow_pressure):
            mdot = (area*self.P*1000)/(np.sqrt(self.R*self.T*1000))*np.sqrt((2*self.k/(self.k-1))*((self.critical_pressure()/self.P)**(2/self.k))*(1-(self.critical_pressure()/self.P)**(1-1/self.k)))
        else:
            mdot = (area*self.P*1000)/(np.sqrt(self.R*self.T*1000))*np.sqrt((2*self.k/(self.k-1))*((backflow_pressure/self.P)**(2/self.k))*(1-(backflow_pressure/self.P)**(1-1/self.k)))

        return mdot

class RocketNozzle(Gas):
    def __init__(self, class_gas):
        self.T=class_gas.T
        self.P=class_gas.P
        self.k=class_gas.k
        self.M=class_gas.M
        self.k=class_gas.k
        self.R=class_gas.R_u/class_gas.M
        self.cp=class_gas.cp
        self.rho = class_gas.P/(class_gas.R*class_gas.T)
        self.normalshock=self.Shock(self)

    def geometry(self, area_start, area_throat, area_end, division=250, color = 'black'):
        A_start = area_start
        A1_star = area_throat
        A_exit  = area_end

        division = 250
        
        r1=int((A_start/A1_star)/(A_start/A1_star+A_exit/A1_star)*division)
        r2=int((A_exit/A1_star)/(A_start/A1_star+A_exit/A1_star)*division)
        
        area_upward = np.linspace((A_start), (A1_star), r1)
        area_downward = np.linspace((A1_star), (A_exit), r2)
        area_total = np.concatenate((area_upward,area_downward))
        
        diameter_total = self.diameter(area_total)
        # plt.style.use('dark_background')
        fig = plt.figure(figsize=(12,6))
        ax = fig.add_subplot(111)
        xs = np.linspace(0,1,r1+r2)
        
        tick_labels=[]     
        for j in np.linspace(0,(r1+r2),11):
            if j==r1+r2:
                tick_labels.append(round(area_total[-1],4))
            else:    
                tick_labels.append(round(area_total[int(j)],4))
        plt.xticks(np.linspace(0,1,11),tick_labels)
        
        plt.plot(xs,diameter_total/2,color=color,linewidth=3)
        plt.plot(xs,-diameter_total/2,color=color,linewidth=3)
        centerline,=plt.plot(xs, 0*xs,linewidth=1,color=color)
        dashes=[30,5,5,5]
        centerline.set_dashes(dashes)
        
        plt.xlabel("Area (m2)")
        plt.ylabel("Radius (m)")
        plt.title("Rocket Nozzle Geometry")
        
        plt.show()
        plt.style.use('default')
        
    def shock(self, exit_pressure, throat_area, exit_area, start_area, plot=True,division = 250):
    
        def shock_finder(A_shock):
            ratio = A_shock/throat_area
            
            M1 = self.ma_finder('downward', ratio)
            P1 = self.pressure(M1, self.P)
            T1 = self.temperature(M1, self.T)
            
            M2 = self.normalshock.Ma2(M1)
            P2 = self.normalshock.P2(M1,P1)
            T2 = self.normalshock.T2(T1, M1)
            
            P02 = self.normalshock.P0_2(self.P, M1)
            
            A2_star = self.normalshock.area_shock_star(throat_area, M1)
            
            ratio2 = exit_area/A2_star
            
            Me = self.ma_finder('upward', ratio2) 
            
            Pe = self.pressure(Me,P02)
            
            target = Pe-exit_pressure
            return target
        if shock_finder(exit_area)>0:
            print("There is no shock wave in the rocket nozzle")
            A_shock = None
        else:
            A_shock=bisection_method( shock_finder,throat_area, exit_area, tolerance = 10e-3,show_iterations=True)
    
        def shock_plot(start_area):
            A_start = start_area
            A1_star = throat_area
            A_exit  = exit_area   
            
            r1=int((A_start/A1_star)/(A_start/A1_star+A_exit/A1_star)*division)
            r2=int((A_exit/A1_star)/(A_start/A1_star+A_exit/A1_star)*division)
            area_upward = np.linspace((start_area), (throat_area), r1)
            area_downward = np.linspace((throat_area), (exit_area), r2)
            area_total = np.concatenate((area_upward,area_downward))

            def find_closest(A, target):
                #A must be sorted
                idx = A.searchsorted(target)
                idx = np.clip(idx, 1, len(A)-1)
                left = A[idx-1]
                right = A[idx]
                idx -= target - left < right - target
                return idx

            idx=find_closest(area_total,A_shock)

            r=self.diameter(A_shock)/2
            plt.style.use('dark_background')

            self.geometry(start_area, throat_area, exit_area,color='white')
            y=np.linspace(r,-r)
            # correction = ((A_shock/throat_area)+(start_area/throat_area))/((exit_area/throat_area)+(start_area/throat_area))
            x=A_shock*np.sin(5000*y)+idx/division
            plt.plot(x,y,color='gold')
            plt.show()
            plt.style.use('default')
    
        if plot==True:
            shock_plot(start_area)
    
        return A_shock
        

""" 

      #############          RELATIONS CLASS         ###########

"""

class relations:
    def change_in_entropy(T2,T1,P2,P1,cp,R):
        return cp*np.log(T2/T1)-R*np.log(P2/P1)





