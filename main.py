
import myplot

# MDA+
import numpy as np

norm_w=False

def dominate(v1, v2):
    for i, j in zip(v1, v2):
        if i < j:
            return False
    return True

def kskyband(data, k):
    do_cnt=np.array([0 for _ in range(len(data))])
    ret=[]
    for i in range(len(data)):
        for j in range(i+1, len(data)):
            if do_cnt[i]>=k:
                break
            if dominate(data[i], data[j]):
                do_cnt[j]+=1
            elif dominate(data[j], data[i]):
                do_cnt[i]+=1
        if do_cnt[i]<k:
            ret.append(i)
    return ret


from itertools import combinations
import heapq


def closeTo(d):
    id_ = int(d)
    idp1 = int(id_ + 1)
    idm1 = int(id_ - 1)
    d1 = abs(d - idm1)
    d2 = abs(d - id_)
    d3 = abs(d - idp1)
    if d1 > d2:
        if d3 > d2:
            return id_
        else:
            return idp1
    else:
        if d3 > d1:
            return idm1
        else:
            return idp1


def rdominate(v1, v2):
    assert (len(v1) == len(v2))
    for i, j in zip(v1, v2):
        if i < j:
            return False
    return True


def vector_len(v):
    ret=0.0
    for i in v:
        ret+=i*i
    return ret**0.5

class cell:

    def __init__(self, b, cur_level, tar_level, K):
        self.bounds = b
        self.dim = int(len(b) / 2)
        self.uHeat = 0.0
        self.rHeat = 0.0
        self.dHeat = 0.0
        self.MaxKMaxUHeat=0.0

        self.cur_l = cur_level
        self.tar_l = tar_level
        self.k = K
        self.vertexes=[]
        self.get_vertexes()
        self.heap = [0.0 for _ in range(K)]
        self.mdap_s_rksyband = []
        self.s_rskyband = []
        self.rkskyband = []
        self.rskyband_lb_MDA = []
        self.cnt=0
        self.children=[]


    def get_vertexes(self):
        sum_lb = 0.0
        for loc_d in range(self.dim):
            sum_lb += self.bounds[loc_d * 2]
        dis = self.bounds[1] - self.bounds[0]
        num_ub = closeTo((1.0 - sum_lb) / dis)
        for cb in combinations(range(self.dim), num_ub):
            v = [0 for _ in range(self.dim)]
            for i in cb:
                v[i] = 1
            one_cb = []
            for i in range(self.dim):
                if v[i] == 1:
                    one_cb.append(self.bounds[2 * i + 1])
                else:
                    one_cb.append(self.bounds[2 * i])
            # TODO delete begin
            # ret = 0.0
            # for v in one_cb:
            #     ret += v * v
            # r = (1.0 / ret) ** 0.5
            # one_cb=[v * r for v in one_cb]
            # TODO delete end
            self.vertexes.append(one_cb)
        if len(self.vertexes)==0:
            print('stop')

    def get_lb_ub(self, p):
        lb = 1e9
        ub = 0.0
        for v in self.vertexes:
            tmp = np.dot(p, v)
            if norm_w:
                tmp/=vector_len(v)
            if tmp > ub:
                ub = tmp
            if tmp < lb:
                lb = tmp
        return lb, ub

    def get_scores(self, ps, P, lb, ub):
        for j in range(len(ps)):
            lb[j], ub[j]=self.get_lb_ub(P[ps[j]])

    def isLeaf(self):
        return self.cur_l >= self.tar_l

    def get_next_children(self, it):
        u = 1.0
        dis = (self.bounds[1] - self.bounds[0]) / 2  # child should be divided by 2
        l = 1.0 - dis * self.dim
        child_p_num = (1 << self.dim)
        while it < child_p_num:
            lv = [0.0 for _ in range(self.dim)]
            iter_cp = it
            for j in range(self.dim):
                if iter_cp % 2 == 0:
                    lv[j] = self.bounds[j * 2]
                else:
                    lv[j] = (self.bounds[j * 2] + self.bounds[j * 2 + 1]) / 2.0
                iter_cp >>= 1
            slv = sum(lv)
            if l < slv < u:
                child_b = [0.0 for _ in range(self.dim * 2)]
                for i in range(self.dim):
                    child_b[i * 2] = lv[i]
                    child_b[i * 2 + 1] = lv[i] + dis
                self.children.append(cell(child_b, self.cur_l + 1, self.tar_l, self.k))
                return self.children[-1], it
            else:
                it += 1
        return None, it

    def MDA_superSet2RKS(self, data):
        theta = self.heap[0]
        tmp = []
        for i in self.s_rskyband:
            ub = i[2]
            if ub > theta:
                tmp.append(i)
        tmp = sorted(tmp, key=lambda item: item[2], reverse=True)
        scores = []
        for i in range(self.k):
            self.rkskyband.append(tmp[i][0])
            self.rskyband_lb_MDA.append(tmp[i][1])
            s = []
            for v in self.vertexes:
                tmp_score=np.dot(v, data[tmp[i][0]])
                if norm_w:
                    tmp_score/=vector_len(v)
                s.append(tmp_score)
            scores.append(s)

        # TODO 这里是对 下文r-dominate count 的BUG的解决方法  begin  2022.2.21
        for i in range(self.k):
            for j in range(i+1, self.k):
                if rdominate(scores[i], scores[j]):
                    self.dHeat+=1
        # TODO 这里是对 下文r-dominate count 的BUG的解决方法 end  2022.2.21
        for i in range(self.k, len(tmp)):
            r_dominate_count = 0
            s = []
            for j in range(len(self.rkskyband)):
                if tmp[i][1] < self.rskyband_lb_MDA[j]:
                    if len(s) == 0:
                        for v in self.vertexes:
                            tmp_score = np.dot(v, data[tmp[i][0]])
                            if norm_w:
                                tmp_score /= vector_len(v)
                            s.append(tmp_score)
                    if rdominate(scores[j], s):
                        r_dominate_count += 1
                        if r_dominate_count >= self.k:
                            break
            if r_dominate_count < self.k:
                # TODO 2022.2.21 这里我们确定原先代码的统计信息错了
                # TODO 这里只统计了其他option对不是前k upper bound option的r-dominate count
                # TODO ANS: 应该在前面补充upper bound前k option他们之间的的r-dominate count
                self.dHeat += r_dominate_count
                self.rkskyband.append(tmp[i][0])
                self.rskyband_lb_MDA.append(tmp[i][1])
                if len(s) == 0:
                    for v in self.vertexes:
                        tmp_score=np.dot(v, data[tmp[i][0]])
                        if norm_w:
                            tmp_score /= vector_len(v)
                        s.append(tmp_score)
                scores.append(s)

    def MDAp_insert(self, parent_sRSK, data, ret):
        lbs = [0.0 for _ in range(len(parent_sRSK))]
        ubs = [0.0 for _ in range(len(parent_sRSK))]
        self.get_scores(parent_sRSK, data, lbs, ubs)
        for lb in lbs:
            theta = self.heap[0]
            if lb > theta:
                heapq.heappush(self.heap, lb)  # heap push
                heapq.heappop(self.heap)
        theta = self.heap[0]
        for j in range(len(parent_sRSK)):
            if ubs[j] > theta:
                self.mdap_s_rksyband.append(parent_sRSK[j])
                if self.isLeaf():
                    self.s_rskyband.append([parent_sRSK[j], lbs[j], ubs[j]])
        if self.cur_l < self.tar_l:
            child_p_num = 1 << self.dim
            j = int(0)
            while j < child_p_num:
                child, j = self.get_next_children(j)
                if child is not None:
                    child.MDAp_insert(self.mdap_s_rksyband, data, ret)
                j += 1
        if self.isLeaf():
            self.MDA_superSet2RKS(data)
            self.uHeat=theta
            # In fact, theta is MaxMin_k since the options with k maximum minimum utility must be in RKS
            # why must be in RKS, since you can't find k option r-dominate them
            # except these k-1 options with maximum minimum utility
            # the other option's minimum score must be lower than k-th option's corresponding score
            self.rHeat=len(self.rkskyband)
            # self.dHeat: updated in self.MDA_superSet2RKS(data)
            ret.append(self)
            if len(ret)%1000==0:
                print(len(ret))

    def __str__(self):
        return self.bounds.__str__()+'['+str(self.uHeat)+','+str(self.rHeat)+','+str(self.dHeat)+']'


def MDAp(root, data, ret):
    idx = [i for i in range(len(data))]
    root.MDAp_insert(idx, data, ret)

from itertools import combinations


def main():
    # read data
    pdt_ta = []
    # with open('./gen_data/cs5_sph1K3d_la.txt', 'r') as f:
    # with open('./hotel.txt', 'r') as f:
    # with open('./NBA-2021-22/NBA2020-21-SF-pure-norm.csv', 'r') as f:
    with open('./NBA-2021-22/NBA2020-21-C-pure-norm.csv', 'r') as f:
        lines = f.readlines()
        for li in lines:
            li = li.strip()
            if li:
                vals = li.split(',') # used for NBA
                # vals = li.split(' ') # used for the others
                vals = [float(val) for val in vals]
                pdt_ta.append(vals)
    # myplot.draw_2d_pdt(pdt_ta, f='./cs5_anti1K3d_la/cs.png')


    print(len(pdt_ta))
    print(len(pdt_ta[0]))
    pdt_ta = np.array(pdt_ta)

    user_ta = []
    with open('./user.txt', 'r') as f:
        lines = f.readlines()
        for li in lines:
            li = li.strip()
            if li:
                vals = li.split(' ')
                vals = [float(val) for val in vals]
                user_ta.append(vals)
    print(len(user_ta))
    print(len(user_ta[0]))
    user_ta = np.array(user_ta)

    # myplot.draw_user(user_ta[:, 0:3], sample=0.02, f='./hotel/TA_user_n_02')
    # exit(0)
    # option dataset to kskyband
    rets=[]
    for attrs in combinations(range(0, 3), 3):
        attrs=[23-5, 24-5, 29-5] # total rebounds TRB:23, assist AST:24, points PTS: 29 this is used for NBA case study
        print('#', attrs)
        k = 10
        d = len(attrs)
        h = 4
        pdt = pdt_ta[:, attrs]
        # myplot.draw_3d_pdt(pdt, f='./NBA-2021-22-out/cs-SF.png')
        # myplot.draw_3d_pdt(pdt, f='./NBA-2021-22-out/cs2-SF.png', view=True)
        kskyid = kskyband(pdt, k)
        pdt_ksky = pdt[kskyid, :]
        print('kskyband size:', len(pdt_ksky))
        b = []
        for _ in range(d):
            b.append(0.0)
            b.append(1.0)
        r = cell(b, 0, h, k)
        ret = []
        MDAp(r, pdt_ksky, ret)
        print("finish algorithm")
        # begin plot and case study
        attrstr = ''
        for a in attrs:
            attrstr += str(a) + '_'
        if norm_w:
            attrstr+='n_'
        # myplot.cs_2d_5cluster(cells=ret, pdt=pdt_ksky)
        # myplot.cs_2d_5cluster(cells=ret, pdt=pdt_ksky, f="./cs5_anti1K3d_la/"+'cs5_K'+str(k)+'H'+str(h), norm=norm_w)
        # myplot.draw_3d_pdt(pdt_ksky, f='./NBA-2021-22-out/sky-SF.png')
        # myplot.draw_3d_pdt(pdt_ksky, f='./NBA-2021-22-out/sky2-SF.png', view=True)
        # myplot.cs_3d_5cluster(cells=ret, pdt=[], f="./NBA-2021-22-out/"+'cs-C_K'+str(k)+'H'+str(h)+'_'+attrstr, norm=norm_w)
        # myplot.cs_3d_5cluster(cells=ret, pdt=pdt_ksky, f="./NBA-2021-22-out/"+'cs6_K'+str(k)+'H'+str(h)+attrstr, norm=norm_w)

        myplot.cs_3d_5cluster(cells=ret, pdt=pdt_ksky, norm=norm_w)

        # myplot.cs1_plot(ret, pdt_ksky, user_ta[:, attrs], attrstr, root=r)
        # myplot.cs2_plot(ret, pdt_ksky, user_ta[:, attrs], attrstr)
        # for c in ret:
        #     print(c.rkskyband, "@", c.bounds)


        rets.append(ret)
    return rets

if __name__ == '__main__':
    lcells=main()
