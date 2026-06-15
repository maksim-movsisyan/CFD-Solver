import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np


def plot_residuals(filename):
    try:
        data = np.loadtxt(filename, skiprows=3).T

        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))
        axes.plot(data[0], data[1], color='royalblue', linewidth=2, label=r'Pressure')
        axes.plot(data[0], data[2], color='orangered', linewidth=2, label=r'$Velocity_x$')
        axes.plot(data[0], data[3], color='green', linewidth=2, label=r'$Velocity_y$')
        axes.plot(data[0], data[4], color='orange', linewidth=2, label=r'$Temperature$')
        axes.set_xlabel('iterations')
        axes.set_ylabel('residuals')
        axes.set_yscale('log', base=10)
        axes.set_title('Convergence history')
        axes.legend()
        axes.grid()
        plt.show()



    except Exception as e:
        print(f"Ошибка: {e}")
        input("Нажмите Enter для выхода...")


if __name__ == "__main__":
    filename='log_file.txt'
    plot_residuals(filename)