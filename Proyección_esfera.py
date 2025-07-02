import numpy as np

class EsferaProyectada:
    def __init__(self, radius=2,
                 n_latitudes_1=8, n_longitudes_1=100,
                 n_latitudes_2=100, n_longitudes_2=8,
                 theta_deg=120, phi_deg=35,
                 include_circulo=True):

        self.radius = radius
        self.n_lat_1 = n_latitudes_1
        self.n_lon_1 = n_longitudes_1
        self.n_lat_2 = n_latitudes_2
        self.n_lon_2 = n_longitudes_2
        self.theta = np.radians(theta_deg)
        self.phi = np.radians(phi_deg)
        self.include_circulo = include_circulo

        # Matriz de rotación
        self.R = self._crear_matriz_rotacion()

    def _crear_matriz_rotacion(self):
        Rz = np.array([
            [np.cos(self.theta), -np.sin(self.theta), 0],
            [np.sin(self.theta),  np.cos(self.theta), 0],
            [0,                   0,                  1]
        ])

        Ry = np.array([
            [ np.cos(self.phi), 0, np.sin(self.phi)],
            [ 0,                1, 0              ],
            [-np.sin(self.phi), 0, np.cos(self.phi)]
        ])

        return Rz @ Ry

    def proyectar(self):
        # Red de latitud y longitud 1
        lat1 = np.linspace(-np.pi/2, np.pi/2, self.n_lat_1)
        lon1 = np.linspace(0, 2*np.pi, self.n_lon_1)
        X1 = self.radius * np.outer(np.cos(lat1), np.cos(lon1))
        Y1 = self.radius * np.outer(np.cos(lat1), np.sin(lon1))
        Z1 = self.radius * np.outer(np.sin(lat1), np.ones_like(lon1))

        # Red de latitud y longitud 2
        lat2 = np.linspace(-np.pi/2, np.pi/2, self.n_lat_2)
        lon2 = np.linspace(0, 2*np.pi, self.n_lon_2)
        X2 = self.radius * np.outer(np.cos(lat2), np.cos(lon2))
        Y2 = self.radius * np.outer(np.cos(lat2), np.sin(lon2))
        Z2 = self.radius * np.outer(np.sin(lat2), np.ones_like(lon2))

        # Aplicar rotación y proyección
        proyecciones = []

        for i in range(self.n_lat_1):
            for j in range(self.n_lon_1):
                x, y, z = self.R @ np.array([X1[i, j], Y1[i, j], Z1[i, j]])
                if y > 0:
                    proyecciones.append([x, z])

        for i in range(self.n_lat_2):
            for j in range(self.n_lon_2):
                x, y, z = self.R @ np.array([X2[i, j], Y2[i, j], Z2[i, j]])
                if y > 0:
                    proyecciones.append([x, z])

        # Círculo de contorno
        if self.include_circulo:
            phi_circle = np.linspace(-np.pi, np.pi, 100)
            for angle in phi_circle:
                proyecciones.append([self.radius * np.cos(angle), self.radius * np.sin(angle)])

        return np.array(proyecciones)


"""import numpy as np
import matplotlib.pyplot as plt

# Parámetros de la esfera
radius = 2
n_latitudes = 10  # Número de líneas de latitud
n_longitudes = 1000  # Número de líneas de longitud

# Crear la esfera
latitudes = np.linspace(-np.pi / 2, np.pi / 2, n_latitudes)
longitudes = np.linspace(0, 2 * np.pi, n_longitudes)

# Coordenadas esféricas a cartesianas
X_1 = radius * np.outer(np.cos(latitudes), np.cos(longitudes))
Y_1 = radius * np.outer(np.cos(latitudes), np.sin(longitudes))
Z_1 = radius * np.outer(np.sin(latitudes), np.ones_like(longitudes))

n_latitudes_2 = 1000  # Número de líneas de latitud
n_longitudes_2 = 10  # Número de líneas de longitud

latitudes_2 = np.linspace(-np.pi / 2, np.pi / 2, n_latitudes_2)
longitudes_2 = np.linspace(0, 2 * np.pi, n_longitudes_2)

X_2 = radius * np.outer(np.cos(latitudes_2), np.cos(longitudes_2))
Y_2 = radius * np.outer(np.cos(latitudes_2), np.sin(longitudes_2))
Z_2 = radius * np.outer(np.sin(latitudes_2), np.ones_like(longitudes_2))

# Dirección aleatoria para la proyección
#theta = np.random.uniform(0, 2 * np.pi)  # Ángulo en el plano XY
#phi = np.random.uniform(0, np.pi)  # Ángulo de altura (zenital)
theta=120/180 * np.pi
phi = 35 /180 * np.pi

Rz = np.array([
    [np.cos(theta), -np.sin(theta), 0],
    [np.sin(theta),  np.cos(theta), 0],
    [0,              0,             1]
])

Ry = np.array([
    [ np.cos(phi), 0, np.sin(phi)],
    [ 0,           1, 0         ],
    [-np.sin(phi), 0, np.cos(phi)]
])

R = Rz @ Ry

view_dir = np.array([np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)])

print(view_dir, R @ np.array([0,0,1]) )

# Proyectar todas las líneas de latitud y longitud


def project_point(x, y, z, view_dir):
    # Proyección estereográfica sobre el plano XY
    dot_product = np.dot(view_dir, np.array([x, y, z]))
    scale_factor = radius / (radius - dot_product)  # Escala de proyección
    return scale_factor * np.array([x, y])


projections = []
for i in range(n_latitudes):
    for j in range(n_longitudes):
        [x, y, z]=R @ np.array([X_1[i,j],Y_1[i,j],Z_1[i,j]]) 
        if y>0:
            projections.append([x,z])


for i in range(n_latitudes_2):
    for j in range(n_longitudes_2):
        [x, y, z]=R @ np.array([X_2[i,j],Y_2[i,j],Z_2[i,j]]) 
        if y>0:
            projections.append([x,z])

projections.append([x,z])
# Convertir a numpy array
N=1000
phi = np.linspace(-np.pi , np.pi , N)

for i in range(N):  
    projections.append([radius * np.cos(phi[i]), radius * np.sin(phi[i])])

projections =np.array(projections)


# Graficar la proyección
plt.figure(figsize=(6, 6))
plt.scatter(projections[:, 0], projections[:, 1], s=1, color='b', alpha=0.7)
plt.title("Proyección de líneas de latitud y longitud")
plt.xlabel("x")
plt.ylabel("y")
#plt.axis('equal')
#plt.grid(True)
plt.show()"""