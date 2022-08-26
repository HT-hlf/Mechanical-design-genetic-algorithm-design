#encoding=UTF-8
import matplotlib.pyplot as plt
import numpy as np

DNA_size = 10
X_BOUND = [0, 5]
cross_rate = 0.8
mutation_rate = 0.003
pop_size = 100
N_GENERATIONS = 200


class GA(object):
    def __init__(self, DNA_size, DNA_bound, cross_rate, mutation_rate, pop_size, n_generations):
        self.DNA_size = DNA_size
        DNA_bound[1] += 1
        self.DNA_bound = DNA_bound
        self.cross_rate = cross_rate
        self.mutate_rate = mutation_rate
        self.pop_size = pop_size
        self.pop = np.random.randint(*DNA_bound, size=(pop_size, self.DNA_size))  # .repeat(pop_size, axis=0)
        # print(self.translateDNA(self. pop))
        self.n_generations = n_generations

    # 模拟函数
    def F(self, x):
        return np.sin(10 * x) * x + np.cos(2 * x) * x

    # DNA翻译
    def translateDNA(self, pop):
        # convert binary DNA to decimal and normalize it to a range(0, 5)
        return pop.dot(2 ** np.arange(self.DNA_size)[::-1]) / float(2 ** self.DNA_size - 1) * X_BOUND[1]

    # 适配
    def get_fitness(self, product):
        return product + 1e-3 - np.min(product)

    # 选择
    def select(self, fitness):
        idx = np.random.choice(np.arange(self.pop_size), size=self.pop_size, replace=True,
                               p=fitness / fitness.sum())
        return self.pop[idx]

    # 交配
    def crossover(self, parent, pop):  # mating process (genes crossover)
        if np.random.rand() < self.cross_rate:
            i_ = np.random.randint(0, self.pop_size, size=1)  # select another individual from pop
            cross_points = np.random.randint(0, 2, size=self.DNA_size).astype(np.bool)  # choose crossover points
            parent[cross_points] = pop[i_, cross_points]  # mating and produce one child
        return parent

    # 突变
    def mutate(self, child):
        for point in range(self.DNA_size):
            if np.random.rand() < self.mutate_rate:
                child[point] = 1 if child[point] == 0 else 0
        return child

    # 进化
    def evolve(self):
        plt.ion()  # something about plotting
        x = np.linspace(*X_BOUND, 200)
        # plt.plot(x, self.F(x))
        for _ in range(self.n_generations):
            F_values = self.F(self.translateDNA(self.pop))  # compute function value by extracting DNA

            # something about plotting
            if 'sca' in locals():
                sca.remove()
            sca = plt.scatter(self.translateDNA(self.pop), F_values, s=200, lw=0, c='red', alpha=0.5);
            plt.pause(5)

            # GA part (evolution)
            fitness = self.get_fitness(F_values)
            print("Most fitted DNA: ", self.pop[np.argmax(fitness), :])
            print(self.translateDNA(self.pop[np.argmax(fitness), :]))
            self.pop = self.select(fitness)
            pop_copy = self.pop.copy()
            for parent in self.pop:
                child = self.crossover(parent, pop_copy)
                child = self.mutate(child)
                parent[:] = child  # parent is replaced by its child

        plt.ioff();
        plt.show()


ga = GA(DNA_size=10,
        DNA_bound=[0, 1],
        cross_rate=0.8,
        mutation_rate=0.05,
        pop_size=100, n_generations=N_GENERATIONS)
ga.evolve()