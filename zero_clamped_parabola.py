import matplotlib.pyplot as plt
import numpy as np

class ZeroClampedParabola:
    """
        Equation: 
            y = max(a(x-b)²+h , 0)
    """
    def __init__(self, x: np.ndarray, a: int, b: int, h: int) -> None:
        self.x = x
        self.a = a
        self.b = b
        self.h = h

        z = a*(x-b)**2 + h
        self.z = np.maximum(0, z)

    def get_consec_point_distances(self) -> np.ndarray:
        x = self.x[:-1]
        z = self.z[:-1]

        x_next = self.x[1:]
        z_next = self.z[1:]

        return np.sqrt((x - x_next)**2 + (z - z_next)**2)

    def get_z_points(self) -> np.ndarray:
        return self.z
    
    def draw_parabola(self, save_file_path: str) -> None:
        plt.plot(self.x, self.z)
        plt.savefig(save_file_path)