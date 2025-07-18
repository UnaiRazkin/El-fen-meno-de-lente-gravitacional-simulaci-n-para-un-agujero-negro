import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from tqdm import tqdm
from PIL import Image
from scipy import optimize
from scipy.special import ellipkinc
from Proyección_esfera import EsferaProyectada
class Num_SchwarzschildRayTracing:
    def __init__(self, M=1.0, L=4.0, E=1.0):
        self.G = 1.0
        self.c = 1.0
        self.M = M
        self.L = L
        self.E = E
        self.rs = 2 * self.G * self.M / self.c**2

    def V(self, r, epsilon):
        return (-0.5 * epsilon +
                epsilon * self.G * self.M / (self.c**2 * r) +
                self.L**2 / (2 * r**2) -
                self.G * self.M * self.L**2 / (self.c**2 * r**3))

    def dV_dr(self, r, epsilon):
        return (-epsilon * self.G * self.M / (self.c**2 * r**2) -
                self.L**2 / r**3 +
                3 * self.G * self.M * self.L**2 / (self.c**2 * r**4))

    def system(self, tau, y, epsilon):
        r, phi, dr_dtau = y

        # Segundo miembro de la ecuación diferencial de segundo orden para r
        d2r_dtau2 = -self.dV_dr(r, epsilon)

        dphi_dtau = self.L / r**2

        return [dr_dtau, dphi_dtau, d2r_dtau2]

    def stop_at_horizon(self, tau, y, epsilon):
        r = y[0]
        return r - self.rs
    stop_at_horizon.terminal = True
    stop_at_horizon.direction = -1

    def solve(self, epsilon, r0=100.0, phi0=0.0, tau_span=(0, 200), max_step=0.1):
        # Resolver el sistema de ecuaciones diferenciales 
        dr0 = -np.sqrt(np.maximum(self.E**2 - 2 * self.V(r0, epsilon), 0))  # Dirección inicial hacia adentro

        y0 = [r0, phi0, dr0]

        sol = solve_ivp(
            fun=lambda tau, y: self.system(tau, y, epsilon),
            t_span=tau_span,
            y0=y0,
            events=lambda tau, y: self.stop_at_horizon(tau, y, epsilon),
            method='RK45',
            max_step=max_step,
            rtol=1e-9,
            atol=1e-9,
            #dense_output=True
            t_eval=np.linspace(*tau_span, 300)
        )

        return sol

class SchwarzschildRayTracer:
    def __init__(self, bh_mass=1.0, Draw_sr_radius=False):

        #self.stars_img = np.array(Image.open("Imágenes/Background_2.png").convert('RGB'))
        #self.stars_img = (np.array(Image.open("Imágenes/Background_3.jpg").convert('RGB')))
        self.stars_img = (np.array(Image.open("Imágenes/bgedit.jpg").convert('RGB')))
        #self.stars_img = np.array(Image.open("Imágenes/Estrellas.jpg").convert('RGB'))

        self.stars_width, self.stars_height = self.stars_img.shape[:2]
        print(self.stars_width, self.stars_height)

        self.image_size_x =  int(self.stars_width/3)
        self.image_size_y =  int(self.stars_height/3)

        self.M = bh_mass
        self.rs = 2 * bh_mass  # Radio de Schwarzschild

        self.r_start = 100 * self.rs
        self.r_max =self.r_start*1.1

        self.b_max_x = 5 * self.rs
        self.b_max_y = self.b_max_x / self.image_size_x * self.image_size_y
        self.b_max = np.sqrt(self.b_max_x**2+self.b_max_y**2)

        self.image = np.zeros((self.image_size_x, self.image_size_y, 3), dtype=np.uint8)  # RGB

        self.draw_sr_radius= Draw_sr_radius
    
    def sphere(self):
        return EsferaProyectada(radius=self.rs).proyectar()
    
    def trace_ray(self, b, theta_view):

        phi= b / self.b_max * np.pi #[0,π]

        def colour(phi, theta_view, alpha):
   
            beta_x = (((phi-alpha) * np.cos(theta_view))/(np.pi)+1)/2 #[0,1]
            beta_y = (((phi-alpha) * np.sin(theta_view))/(np.pi)+1)/2 #[0,1]

            x_img = int(beta_x * self.stars_width) % self.stars_width #[0, 1, ..., self.stars_width]
            y_img = int(beta_y * self.stars_height) % self.stars_height

            return self.stars_img[x_img, y_img]
        
        def draw_sphere(x,y,z):
            #No es así es solamente una prueba

            X, Z, Y= self.sphere()

            for i in range(len(X)):
                for j in range(len(Y)):
                    for k in range(len(Z)):
                        if (z<0 and Z[k]<0):
                            if self.closeness(X[i], Y[j], x, y):
                                return [0,0,255]
                        elif (z>0 and Z[k]>0):
                            if self.closeness(X[i], Y[j], x, y):
                                return [0,0,255]
            print("bien")
            return [0,0,0]
        rm=b
        def rmin(R, b):
            return 1 / b**2 - 1 / R**2 + self.rs / R**3

        r_min = self.rs * 3/2
        r_max = self.r_max

        fa = rmin(r_min, b)
        fb = rmin(r_max, b)

        if fa * fb > 0:
            if b>3/2 * self.rs:
                fa = rmin(self.rs, b)
            else:
                #print("bien")
                tracer_massless = Num_SchwarzschildRayTracing(M=self.M, L=b, E=1)
                sol_light = tracer_massless.solve(epsilon=0, r0=self.r_start, max_step=1, tau_span=(0, 200))
                x = sol_light.y[0] * np.sin(sol_light.y[1]) * np.cos(theta_view)
                y = sol_light.y[0] * np.sin(sol_light.y[1]) * np.sin(theta_view)
                z = sol_light.y[0] * np.cos(sol_light.y[2]) 
                return draw_sphere(x[-1],y[-1],z[-1])
        else:
            sol = optimize.root_scalar(rmin, bracket=[r_min, r_max], args=(b,), method='brentq')
            rm = sol.root
        
        s=np.sqrt((rm-self.rs)*(rm+3*self.rs))
        m=(s-rm+3*self.rs)/2/s
        arg=np.sqrt(2*s/(3*rm-3*self.rs+s))
        arg = np.clip(arg, -1.0, 1.0)
        varphi=np.arcsin(arg)
        
        alpha=np.abs((4*np.sqrt(rm/s)*float(ellipkinc(varphi, m)) -np.pi)) % (np.pi)   #[0, π]
        #alpha=0
        if np.isnan(alpha):
            return [0,0,0]
        else:
            return colour(phi, theta_view, alpha)

    def closeness(self, X, Y, x, y):

        b_max_x = self.b_max_x
        b_max_y = self.b_max_y

        A= ((X - 2 * b_max_x / (self.image_size_x / 2) < x) and (x < X + 2 * b_max_x / (self.image_size_x / 2) ))
        B= ((Y - 2 * b_max_y / (self.image_size_y / 2) < y) and (y < Y + 2 * b_max_y / (self.image_size_y / 2) ))

        return (A and B)
    
    def render(self):
        
        d=np.sqrt(self.image_size_x**2/self.image_size_y**2+1)

        with tqdm(total=self.image_size_y, desc="Renderizando", unit="px") as pbar:

            for j in range(self.image_size_y):
                for i in range(self.image_size_x):

                    b_max_x = self.b_max_x
                    b_max_y = self.b_max_y

                    x = (i - self.image_size_x / 2) / (self.image_size_x / 2) * b_max_x
                    y = (j - self.image_size_y / 2) / (self.image_size_y / 2) * b_max_y
                    
                    b = np.sqrt(x**2 + y**2)

                    theta_view = np.arctan2(y,x)

                    #b_min= self.M / (3* np.sqrt(3))

                    if self.draw_sr_radius:

                        X_rs = self.rs * np.cos(theta_view)
                        Y_rs = self.rs * np.sin(theta_view)
                        
                        X_rphoton = 3/2*self.rs * np.cos(theta_view)
                        Y_rphoton = 3/2*self.rs * np.sin(theta_view)

                        if self.closeness(X_rs , Y_rs , x, y):
                            self.image[i, j] = [255,0,0]

                        elif self.closeness(X_rphoton , Y_rphoton , x, y):
                            self.image[i, j] = [0,255,0]
                        else:
                            color = self.trace_ray(b, theta_view)
                            self.image[i, j] = color

                    #if np.all(self.image[i, j] == [0,0,0]):
                    else:
                        color = self.trace_ray(b, theta_view)
                        self.image[i, j] = color
                pbar.update(1)
        

        #plt.imshow(self.stars_img, cmap='gray', origin='upper')
        plt.imshow(self.image, cmap='gray', origin='upper')
        #plt.title("Sombra del agujero negro (modelo simple)")
        plt.axis('off')
        #plt.show()
        plt.savefig("../Figuras/Sombra del agujero negro.png")

rt = SchwarzschildRayTracer(bh_mass=1.0, Draw_sr_radius=True)
rt.render() # máximo parámetro de impacto

