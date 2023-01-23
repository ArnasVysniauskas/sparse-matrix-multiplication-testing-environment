import csv
from dataclasses import dataclass, field
import math
from typing import DefaultDict, NewType

import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

DATA_FILE = "results_gustavson_blas_v1.csv"

AlgorithmResult = NewType("AlgorithmResult", dict[tuple[int, float], list[int]])
AlgorithmAverage = NewType("AlgorithmAverage", dict[tuple[int, float], float])
AlgorithmDeviation = NewType("AlgorithmDeviation", dict[tuple[int, float], float])
AlgorithmData = (AlgorithmResult | AlgorithmAverage | AlgorithmDeviation)
AlgorithmName = NewType("AlgorithmName", str)

def sign(value: float) -> int:
    if value >= 0: return 1
    return -1

def get_data_from_file(file_name: str) -> tuple[dict[AlgorithmName, AlgorithmResult], list[int], list[float]]:
    algorithm_results: dict[AlgorithmName, AlgorithmResult] = {}
    dimensions: set[int] | list[int] = set() # set to discard duplicates
    sparsities: set[float] | list[float] = set() # set to discard duplicates
    with open(file_name, newline='') as csvfile:
        data_reader = csv.reader(csvfile, delimiter=',')

        algorithm_names: list[AlgorithmName] = [name.strip() for name in next(data_reader)[2:]]
        for name in algorithm_names:
            algorithm_results[name] = DefaultDict(list)

        for line in data_reader:
            dimension = int(line[0])
            sparsity = float(line[1])
            values = [int(val) for val in line[2:]]

            dimensions.add(dimension)
            sparsities.add(sparsity)
            for idx, name in enumerate(algorithm_names):
                key: tuple[int, float] = tuple([dimension, sparsity])
                algorithm_results[name][key].append(values[idx])

    dimensions = sorted(dimensions)
    sparsities = sorted(sparsities)
    return algorithm_results, dimensions, sparsities

def get_algorithm_averages(
    algorithm_results: dict[AlgorithmName, AlgorithmResult], sample: int = -1
) -> dict[AlgorithmName, AlgorithmAverage]:
    algorithm_averages: dict[AlgorithmName, AlgorithmAverage] = {}
    for name in algorithm_results.keys():
        algorithm_averages[name] = {}
        for key, all_values in algorithm_results[name].items():
            algorithm_averages[name][key] = sum(all_values[:sample]) / len(all_values[:sample])
    return algorithm_averages

def get_algorithm_deviations(
    algorithm_results: dict[AlgorithmName, AlgorithmResult],
    algorithm_averages: dict[AlgorithmName, AlgorithmAverage],
    sample: int = -1
) -> dict[AlgorithmName, AlgorithmDeviation]:
    algorithm_deviations: dict[AlgorithmName, AlgorithmDeviation] = {}
    for name in algorithm_results.keys():
        algorithm_deviations[name] = {}
        for key, all_values in algorithm_results[name].items():
            current_average = algorithm_averages[name][key]
            normalised_values = [(current_average - value)**2 for value in all_values[:sample]]
            algorithm_deviations[name][key] = math.sqrt(sum(normalised_values)) / len(all_values[:sample])
    return algorithm_deviations

def get_planes_intersection(
    algorithm_averages: dict[AlgorithmName, AlgorithmAverage], 
    dimensions: list[int], sparsities: list[float],
    algorithm_1: str, algorithm_2: str
) -> tuple[list[int], list[float], list[float]]:
    algorithm_1_results = algorithm_averages[algorithm_1]
    algorithm_2_results = algorithm_averages[algorithm_2]
    
    dimension_intersections = []
    sparsity_intersections = []
    value_intersections = []
    for d in dimensions:
        algorithm_1_values_for_d = [algorithm_1_results[tuple([d, s])] for s in sparsities]
        algorithm_2_values_for_d = [algorithm_2_results[tuple([d, s])] for s in sparsities]
        diff = [v1 - v2 for v1, v2 in zip(algorithm_1_values_for_d, algorithm_2_values_for_d)]
        
        current_sign = sign(diff[0])
        for idx, v in enumerate(diff):
            new_sign = sign(v)
            if new_sign != current_sign:
                current_sign = new_sign
                sparsity = (sparsities[idx - 1] + sparsities[idx]) / 2
                sparsity_intersections.append(sparsity)
                dimension_intersections.append(d)
                value_intersections.append(algorithm_1_values_for_d[idx] - v / 2)
                break
    return dimension_intersections, sparsity_intersections, value_intersections

def get_mesh(
    algorithm_name: AlgorithmName,
    dimensions_mesh: np.ndarray,
    sparsities_mesh: np.ndarray,
    values: dict[AlgorithmName, AlgorithmData]
) -> np.ndarray:
    shape = dimensions_mesh.shape
    result_mesh = np.zeros(shape)
    for idx_c in range(shape[0]):
        for idx_r in range(shape[1]):
            result_mesh[idx_c][idx_r] = values[algorithm_name][
                tuple([dimensions_mesh[idx_c][idx_r], sparsities_mesh[idx_c][idx_r]])
            ]
    return result_mesh

@dataclass
class Plot:
    plot_type: str
    pos: int
    zlabel: str
    title: str
    color: str

    xdata: any
    ydata: any
    zdata: any

@dataclass
class Plotter:

    width: int
    height: int
    view_angle = [45, -120, 0]
    __plots: list[Plot] = field(default_factory=list) 

    def add(self, plot_type: str, column: int, row: int, zlabel: str, title: str, xdata, ydata, zdata = None):
        self.__plots.append(Plot(
            plot_type=plot_type,
            pos = column + row * self.width + 1,
            zlabel=zlabel,
            title=title,
            color="red",
            xdata=xdata,
            ydata=ydata,
            zdata=zdata
        ))

    def plot_all(self):
        fig = plt.figure(figsize=plt.figaspect(1/self.width))
        self.__plots.sort(key = lambda x: x.pos)
        current_pos = -1

        for plot in self.__plots:
            if plot.pos != current_pos:
                current_pos = plot.pos
                ax = fig.add_subplot(self.height, self.width, current_pos, projection='3d')
                ax.set_ylabel("sparsity")
                ax.set_xlabel("dimension")
                ax.set_zlabel(plot.zlabel)
                ax.view_init(self.view_angle[0], self.view_angle[1], self.view_angle[2])
                ax.set_title(plot.title)

                ax.xaxis.label.set_fontsize(10)
                ax.yaxis.label.set_fontsize(10)
                ax.title.set_fontsize(10)

            match plot.plot_type:
                case "plot":
                    ax.plot(
                        plot.xdata, plot.ydata, plot.zdata,
                        color="black", marker ="."
                    )
                case "wireframe":
                    ax.plot_wireframe(
                        plot.xdata, plot.ydata, plot.zdata,
                        color=plot.color, linewidth=.7
                    )
                case "surface":
                    surf = ax.plot_surface(
                        plot.xdata, plot.ydata, plot.zdata,
                        cmap=cm.coolwarm
                    )
                    fig.colorbar(surf, shrink=0.5, aspect=5)
        plt.show()

def main():
    algorithm_results, dimensions, sparsities = get_data_from_file(DATA_FILE)
    no_algorithms = len(algorithm_results.keys())
    algorithm_averages = get_algorithm_averages(algorithm_results)
    algorithm_deviations = get_algorithm_deviations(
        algorithm_results, algorithm_averages
    )
    
    (intersection_dimensions, intersection_sparsities, intersection_values
    ) = get_planes_intersection(
        algorithm_averages, dimensions, sparsities,
        list(algorithm_results.keys())[0], list(algorithm_results.keys())[1]
    )

    # Convert to format acceptable by plotter
    dimensions, sparsities = np.meshgrid(dimensions, sparsities)

    plotter = Plotter(3, no_algorithms)
    for idx, algorithm_name in enumerate(algorithm_results.keys()):
        plotter.add(
            "plot", 0, idx, "computation_time (a.u.)",
            f"computation time of {algorithm_name}",
            intersection_dimensions, intersection_sparsities, intersection_values
        )
        plotter.add(
            "wireframe", 0, idx, "computation_time (a.u.)",
            f"computation time of {algorithm_name}",
            dimensions, sparsities,
            get_mesh(algorithm_name, dimensions, sparsities, algorithm_averages)
        )
        plotter.add(
            "surface", 1, idx, "computation_time std. deviation (a.u.)",
            f"time std. deviation of {algorithm_name}",
            dimensions, sparsities,
            get_mesh(algorithm_name, dimensions, sparsities, algorithm_deviations)
        )
        plotter.add(
            "surface", 2, idx, "ratio of deviation and computation time",
            f"deviation ratio for {algorithm_name}",
            dimensions, sparsities,
            get_mesh(algorithm_name, dimensions, sparsities, algorithm_deviations) / get_mesh(algorithm_name, dimensions, sparsities, algorithm_averages)
        )
    plotter.plot_all()

def deviation_analysis():
    algorithm_results, dimensions, sparsities = get_data_from_file(DATA_FILE)
    no_algorithms = len(algorithm_results.keys())
    algorithms = list(algorithm_results.keys())

    algorithm_averages = []
    algorithm_deviations = []

    algorithm_averages.append(get_algorithm_averages(algorithm_results, 2))
    algorithm_averages.append(get_algorithm_averages(algorithm_results, 3))
    algorithm_averages.append(get_algorithm_averages(algorithm_results, 4))
    algorithm_averages.append(get_algorithm_averages(algorithm_results, 5))
    algorithm_averages.append(get_algorithm_averages(algorithm_results, 6))
    algorithm_averages.append(get_algorithm_averages(algorithm_results, 7))
    algorithm_averages.append(get_algorithm_averages(algorithm_results, 8))
    algorithm_averages.append(get_algorithm_averages(algorithm_results, 9))
    algorithm_averages.append(get_algorithm_averages(algorithm_results, 10))

    algorithm_deviations.append(get_algorithm_deviations(algorithm_results, algorithm_averages[0], 2))
    algorithm_deviations.append(get_algorithm_deviations(algorithm_results, algorithm_averages[1], 3))
    algorithm_deviations.append(get_algorithm_deviations(algorithm_results, algorithm_averages[2], 4))
    algorithm_deviations.append(get_algorithm_deviations(algorithm_results, algorithm_averages[3], 5))
    algorithm_deviations.append(get_algorithm_deviations(algorithm_results, algorithm_averages[4], 6))
    algorithm_deviations.append(get_algorithm_deviations(algorithm_results, algorithm_averages[5], 7))
    algorithm_deviations.append(get_algorithm_deviations(algorithm_results, algorithm_averages[6], 8))
    algorithm_deviations.append(get_algorithm_deviations(algorithm_results, algorithm_averages[7], 9))
    algorithm_deviations.append(get_algorithm_deviations(algorithm_results, algorithm_averages[8], 10))

    dimensions, sparsities = np.meshgrid(dimensions, sparsities)
    plotter = Plotter(3, 3)
    for idx, deviations in enumerate(algorithm_deviations):
        plotter.add(
            "surface", idx, 0, "computation_time std. deviation (a.u.)",
            f"time std. deviation of {algorithms[0]}, sample = {idx + 2}",
            dimensions, sparsities,
            get_mesh(algorithms[0], dimensions, sparsities, deviations)
        )
    plotter.plot_all()


if __name__ == "__main__":
    deviation_analysis()