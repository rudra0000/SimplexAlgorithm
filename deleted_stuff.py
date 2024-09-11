    # def find_q(self):  # find column in non-basis to swap
        # min_ind = 0
        # min_val = None
        # for i in range(1, self.d + 1):
        #     if min_val is None:
        #         min_val = self.Matrix[self.c + 1][i]
        #         min_ind = i
        #     if self.Matrix[self.c + 1][i] < min_val and (
        #         self.is_basis.get(i) == None or self.is_basis.get(i) == False
        #     ):
        #         min_val = self.Matrix[self.c + 1][i]
        #         min_ind = i
        # if min_ind == 0 or min_val >= 0:
        #     return [False, -1]
        # return [True, min_ind]