import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import math as m
import matplotlib.animation as animation
from celluloid import Camera
import ffmpeg

#свойства воды
lam = 0.6
Cp = 4200
ro = 1000

L = 1                           #длина рассчитываемого участка, м. При 2,8 Число Куранта примерно 1
nx = 200                         #число разбиений по длине
nt = 800                       #число разбиений по времени
dx=L/nx                         #размер одной ячейки, м
dtmax = 0.5*ro*Cp*dx*dx/lam     #макс шаг инт по времени, с
V = 50                          #скорость воды, м/с
dt = dtmax*0.000001                             #шаг инт. по времени
RoSkach = 5000

#задание матриц линейных разностных операторов
gam = round(V*dt/dx,4)
B1 = np.eye(nt)
B0 = np.zeros((nt,nt))

for i in range(nt-1):
    for j in range(nx-1):
        if i == j:
            B0[i][j] = 1 - gam
        B0[1][1] = 1
        if i == j+1:
            B0[i][j] = gam

B1obr = np.linalg.inv(B1) #Обратная матрица В1

#Оператор перехода (шага)
def S(yn):
    s = B1obr.dot(B0)
    return s.dot(yn)

RaspRo = np.zeros((nt,nx-1))      #введение матрицы плотностей
RaspRo[0][0] = RoSkach          #величина толчка давления, нач. усл.

def RhoPlusDt(RhoI,RhoIMinusOne,V): #связь плотностей через время dt
    return RhoI-(dt*V/dx)*(RhoI-RhoIMinusOne)

#for j in range(1,nt): #заполнение матрицы плотностей БЕЗ МАТРИЦ
#    for i in range(1,nx-1):
#            RaspRo[i][j] = round(RhoPlusDt(RaspRo[i][j-1], RaspRo[i-1][j-1],RaspV[i+1][j]),2)

for i in range(1,nt):
    RaspRo[i]=S(RaspRo[i-1])

def RaspRoFunc(delay):
    for j in range(0,nx):
        y[j] = RaspRo[j][delay]
    return y

# Функция, вызываемая для каждого кадра
def main_func(frame, line, x):
    '''
    frame - параметр, который изменяется от кадра к кадру.
    line - кривая, для которой изменяются данные.
    x - список точек по оси X, для которых рассчитывается функция Гаусса.
    sigma - отвечает за ширину функции Гаусса.
    '''
    y = RaspRoFunc(frame)
    line.set_ydata(y)
    return [line]

if __name__ == '__main__':
    # Параметры отображаемой функции
    maxSize = 1000
    # Диапазон точек для расчета графика функции
    x = np.arange(maxSize)/nx
    y = np.zeros(maxSize)
    # Создание окна для графика
    fig, ax = plt.subplots()
    # Установка отображаемых интервалов по осям
    ax.set_xlim(0, L)
    ax.set_ylim(0, 1000)
    # Создание линии, которую будем анимировать
    line, = ax.plot(x, y)
    # !!! Параметр, который будет меняться от кадра к кадру
    frames = np.arange(0, nt)
    # !!! Задержка между кадрами в мс
    interval = 20
    # !!! Использовать ли буферизацию для устранения мерцания
    blit = True
    # !!! Будет ли анимация циклической
    repeat = True
    # !!! Создание анимации
    animation = FuncAnimation(
            fig,
            func=main_func,
            frames=frames,
            fargs=(line, x),
            interval=interval,
            blit=blit,
            repeat=repeat)
    plt.show()

print('Число Куранта = ',dt*V/dx)
print('dt = ', dt)
