import matplotlib.pyplot as plt 
import numpy as np

class Plotter:
    def __init__(self, epsilon):
        self.epsilon = epsilon
        self.counter = 0
        plt.ion() ## Note this correction
        fig=plt.figure()
        plt.axis([0,20,0,10])
        self.uxlist = []
        self.uylist = []
        self.oxlist = []
        self.oylist = []
        plt.axhline(y=(1-self.epsilon), color='r', linewidth=4)
        plt.draw()
        plt.pause(0.001)

    def draw(self, y1, y2):
        self.counter += 1
        self.__draw(self.counter, y1, under=True)
        self.__draw(self.counter, y2, under=False)

    def __draw(self, x, y, under=True):
        if under:
            self.uxlist.append(x)
            self.uylist.append(y)
            
            xlist = self.uxlist
            ylist = self.uylist
            color = 'bo'
        else:
            self.oxlist.append(x)
            self.oylist.append(y)

            xlist = self.oxlist
            ylist = self.oylist
            color = 'ro'

        plt.plot(xlist, ylist, color, xlist, ylist, 'k')
        plt.axhline(y=(1-self.epsilon), color='r', linewidth=4)
        plt.draw()
        plt.pause(0.001)
