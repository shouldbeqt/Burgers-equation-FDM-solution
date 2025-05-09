import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import io
from PIL import Image


viscosity = [1, 0.1, 0.01, 0.001]
dx =  10
dt =  100
X = 1
T = 1
x = np.linspace (0, X, dx)
t = np.linspace (0, T, dt)
delta_x = X/(dx-1)
delta_t = T/(dt-1)

class Main_attribute ():
    """
    Class initializate for main attribute which will be using in all mesh
    """
    def __init__ (self):
        self.dx =  10
        self.dt =  1000

        self.X = 1
        self.T = 1

        self.x = np.linspace (0, self.X, self.dx)
        self.t = np.linspace (0, self.T, self.dt)

        self.delta_x = self.X/(self.dx-1)
        self.delta_t = self.T/(self.dt-1)


    def condition (self, matrix):
        """
        U(0, x) = sin (pi*x)
        U(t, 0) = 0
        U(t, L) = 0
        """
        matrix [0, :] = np.sin(np.pi * self.x)
        matrix [:, 0] = 0
        matrix [:, -1] = 0
        return matrix
    
    def plot_visualization_1d(self, matrix):
        count = 0
        data = []
        for i in matrix:
            count += 1
            if count % (self.dt/10) == 0:
                time_val = count / self.dt
                for x_val, y_val in zip(self.x, i):
                    data.append({"x": x_val, "y": y_val, "time": time_val})
        df = pd.DataFrame(data)
        plt.figure(figsize=(10, 6))
        sns.lineplot(data=df, x="x", y="y", hue="time", palette="viridis")
        plt.title(self.type)
        plt.xlabel("Координата x")
        plt.ylabel("Скорость u(x,t)")
        plt.legend(title="время", loc='upper right', bbox_to_anchor=(1, 1))
        plt.tight_layout()
        plt.savefig(f"1d_mesh_{self.type}_viscosity_{self.viscosity}.png", format='png') 
        plt.close()
        

    def contour_visualization (self, matrix):
        plt.contourf (self.t, self.x, matrix.T, cmap = 'plasma')
        cbar = plt.colorbar(aspect=10, shrink=0.8)
        cbar.set_label('Value')
        plt.title(self.type)
        plt.ylabel("Координата x")
        plt.xlabel("Время t")
        plt.savefig(f"2d_mesh_{self.type}_viscosity_{self.viscosity}.png", format='png') 
        plt.close()



        
    def CFL_requirments (self):
        if self.delta_t/self.delta_x**2 >= 0.5:
            raise ValueError ("CFL requirments frustrate")
        else:
            print ("Mesh sustainable")

class Straight_mesh (Main_attribute):
    def __init__ (self, viscosity):
        super().__init__()
        self.type = "Явная схема против потока"
        self.main = np.zeros ( (self.dt, self.dx) )
        self.main = self.condition (self.main)
        self.viscosity = viscosity
        self.main = self.explicit_mesh (self.main)
        self.plot_visualization_1d (self.main)
        self.contour_visualization (self.main)



    def explicit_mesh (self, matrix):
        self.CFL_requirments ()
        for i in range (self.dt-1):
            for j in range (1, self.dx-1):
                diffusion = self.viscosity * (matrix[i, j+1] - 2*matrix[i, j] + matrix[i, j-1]) / self.delta_x**2
                convection = matrix[i, j] * (matrix[i, j] - matrix[i, j-1]) / self.delta_x
                matrix[i+1, j] = matrix[i, j] + self.delta_t * (diffusion - convection)

        return matrix
    
    

    
class Laks_mesh (Main_attribute):
    def __init__ (self, viscosity):
        super().__init__()
        self.type = "Схема Лакса"
        self.main = np.zeros ( (self.dt, self.dx) )
        self.main = self.condition (self.main)
        self.viscosity = viscosity
        self.main = self.laks_mesh (self.main)
        
        self.plot_visualization_1d (self.main)
        self.contour_visualization (self.main)


    def laks_mesh (self, matrix):
        self.CFL_requirments ()
        for i in range (self.dt-1):
            for j in range (1, self.dx-1):
                laks_part = 0.5 * (matrix[i, j+1] + matrix [i, j-1])
                diffusion = self.viscosity * (matrix[i, j+1] - 2*matrix[i, j] + matrix[i, j-1]) / self.delta_x**2
                convection = matrix[i,j] * (matrix[i, j+1] - matrix[i, j-1]) / (2*self.delta_x)
                
                matrix[i+1, j] = self.delta_t * (diffusion - convection) + laks_part
                
            
        return matrix

class Makkormak_mesh (Main_attribute):
    def __init__(self, viscosity):
        super().__init__() 
        self.type = "Схема Маккормака"
        self.main = np.zeros ((self.dt, self.dx))
        self.main = self.condition (self.main)
        self.viscosity = viscosity
        self.main = self.makkormak_mesh (self.main)
        self.plot_visualization_1d (self.main)
        self.contour_visualization (self.main)

    def makkormak_mesh(self, matrix):
        self.CFL_requirments ()

        predictor = matrix.copy()
        corrector = matrix.copy()

        for i in range(self.dt - 1):

            for j in range(1, self.dx - 1):
                conv = matrix[i, j] * (matrix[i, j+1] - matrix[i, j]) / self.delta_x
                diff = self.viscosity * (matrix[i, j+1] - 2*matrix[i, j] + matrix[i, j-1]) / self.delta_x**2
                predictor[i+1, j] = matrix[i, j] + self.delta_t * (diff - conv)

            for j in range(1, self.dx - 1):
                conv_corr = predictor[i+1, j] * (predictor[i+1, j] - predictor[i+1, j-1]) / self.delta_x
                diff_corr = self.viscosity * (predictor[i+1, j+1] - 2*predictor[i+1, j] + predictor[i+1, j-1]) / self.delta_x**2
                corrector[i+1, j] = matrix[i, j] + self.delta_t * (diff_corr - conv_corr)

            for j in range(1, self.dx - 1):
                matrix[i+1, j] = 0.5 * (predictor[i+1, j] + corrector[i+1, j])

        return matrix
                
class Lacks_Vendrof (Main_attribute):
    def __init__(self, viscosity):
        super().__init__() 
        self.type = "Схема Лакса - Вендроффа"
        self.main = np.zeros ((self.dt, self.dx))
        self.main = self.condition (self.main)
        self.viscosity = viscosity
        self.main = self.lacks_vendrof_mesh (self.main)
        self.plot_visualization_1d (self.main)
        self.contour_visualization (self.main)

    def lacks_vendrof_mesh (self, matrix):
        self.CFL_requirments ()
        for i in range (self.dt-1):
            predictor = np.zeros (self.dx)
            for j in range (1, self.dx-1):
                conv = matrix[i, j] * (matrix[i, j] - matrix[i, j-1]) / (self.delta_x)
                diff = self.viscosity * (matrix[i, j + 1] - 2 * matrix[i, j] + matrix[i, j - 1]) / (self.delta_x ** 2)
                predictor[j] = 0.5 * (matrix[i, j+1]+ matrix[i,j]) + (-conv + diff) * (self.delta_t / 2) 

            for j in range (1, self.dx-1):
                conv_predictor = predictor[j] * (predictor[j + 1] - predictor[j - 1]) / (self.delta_x)
                diff_predictor = self.viscosity * (predictor[j + 1] - 2 * predictor[j] + predictor[j - 1]) / (self.delta_x ** 2)
                matrix[i + 1, j] = matrix[i, j] + (-conv_predictor + diff_predictor) * self.delta_t
                
        return matrix

def create_graph (array_1, array_2, array_3, array_4, name):
        plt.clf ()
        plt.plot (x, array_1, label = 'Вязкость: 1')
        plt.plot (x, array_2,label = 'Вязкость: 0.1')
        plt.plot (x, array_3,label = 'Вязкость: 0.01')
        plt.plot (x, array_4, label = 'Вязкость: 0.001')
        plt.title (name)
        plt.legend (title="Вязкость", loc='upper right', bbox_to_anchor=(1, 1))
        plt.xlabel("Координата x")
        plt.ylabel("Скорость u(x,t)")
        plt.xlim (0, 1)
        plt.ylim (0, 1)

def create_gif (main_matrix, name):
    frames = []
    for i in range (dt):
            if i % (dt//100) == 0:
                create_graph (main_matrix[0][i], main_matrix[1][i], main_matrix[2][i], main_matrix[3][i], name)
                buf = io.BytesIO()
                plt.savefig(buf, format='png')
                buf.seek(0)
                frames.append(Image.open(buf))
    frames[0].save(f'{name}.gif',
              save_all=True,
              append_images=frames[1:],
              duration=10,
              loop=0)


    
if __name__ == "__main__":
    names = ["Явная схема против потока", "Схема Лакса", "Схема Маккормака", "Схема Лакса-Вендроффа"]
    mesh_types = [Straight_mesh, Laks_mesh, Makkormak_mesh, Lacks_Vendrof]
    for mesh_index, mesh_class in enumerate(mesh_types):
        for_gif = []
        for temp_visc in viscosity:
            temp_mesh = mesh_class(temp_visc)
            print (temp_mesh.type)
            for_gif.append(temp_mesh.main)
        
        print(f"Создаем GIF для {names[mesh_index]} с {len(for_gif)} наборами данных")
        create_gif(for_gif, names[mesh_index])
    
