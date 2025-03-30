import numpy as np
from zero_clamped_parabola import ZeroClampedParabola

class SkinModelConstParams:
    def __init__(self, xlimits: list, ylimits: list, zlimits: list, num_x_ticks: int, num_y_ticks: int, num_z_ticks: int, tx_rx_distance: int, mean_penetration: int) -> None:
        self.x = np.linspace(xlimits[0], xlimits[1], num_x_ticks)
        self.y = np.linspace(ylimits[0], ylimits[1], num_y_ticks)
        self.z = np.linspace(zlimits[0], zlimits[1], num_z_ticks)
        self.meshY, self.meshZ, self.meshX = np.meshgrid(self.y.astype(np.single), self.z.astype(np.single), self.x.astype(np.single))

        voxel_len_y = self.y[1] - self.y[0]
        voxel_len_z = self.z[1] - self.z[0]
        self.voxel_yz_area = voxel_len_y * voxel_len_z

        a = -(4 * mean_penetration) / (tx_rx_distance ** 2)
        parabola = ZeroClampedParabola(self.x, a, 0, mean_penetration)
        self.mean_light_pathway = parabola
        self.mean_light_path_dist = parabola.get_consec_point_distances()
        self.mean_light_pathway_valid_indices = (self.x > -tx_rx_distance / 2) & (self.x < tx_rx_distance / 2)

        self.ln10 = np.log(10)
    
    def get_mean_light_path_distances(self) -> np.ndarray:
        return self.mean_light_path_dist
    
    def get_mean_light_path_z_vals(self) -> np.ndarray:
        return self.mean_light_pathway.get_z_points()
    
    def is_valid_mean_light_path_index(self, index: int) -> bool:
        return self.mean_light_pathway_valid_indices[index]
    
    def get_meshX(self) -> np.ndarray:
        return self.meshX
    
    def get_meshY(self) -> np.ndarray:
        return self.meshY
    
    def get_meshZ(self) -> np.ndarray:
        return self.meshZ
    
    def get_x(self) -> np.ndarray:
        return self.x
    
    def get_y(self) -> np.ndarray:
        return self.y
    
    def get_z(self) -> np.ndarray:
        return self.z
    
    def get_voxel_yz_area(self) -> int:
        return self.voxel_yz_area

    def get_ln10(self) -> float:
        return self.ln10
    

    
