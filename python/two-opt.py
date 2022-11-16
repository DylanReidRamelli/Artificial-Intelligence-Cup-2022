class TwoOpt:
    def __init__(self, tour):
        self.tour = tour
        self.tsp_solution = []
        self.best_gain = 0
        
    def compute_gain(i,j):
        # Compute the gain obtained frm performing 
        return 0
        
    def run(self):
        while self.best_gain != 0:
            self.best_gain = 0
            for i in range(1,len(self.tour)):
                for j in range (0,len(self.tour)):
                    gain = self.compute_gain(i,j)
                    if gain < self.best_gain:
                        self.best_gain = gain
