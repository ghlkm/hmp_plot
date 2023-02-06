
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.colors as colors
from matplotlib import cm
from scipy.spatial import distance
import dictances
import pandas as pd
import matplotlib


import matplotlib.pyplot as plt
# this should be set for all figures
# since SIGMOD requests all figures should be True-Type or Type-1 font embedding
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def bhattacharyya_coefficient(v1, v2):
    v1s=sum(v1)
    v2s=sum(v2)
    v1n=[i/v1s for i in v1]
    v2n=[i/v2s for i in v2]
    ret=0.0
    for p, q in zip(v1n, v2n):
        ret+=(p*q)**0.5
    return ret

def my3Dheatmap(cells, heats, fn='1'):
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    # fig = plt.figure()
    fig = plt.figure(figsize=(16, 16))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_box_aspect([1, 1, 1])
    ax.view_init(azim=45, elev=90*np.arccos(1.0/(3.0**0.5))/(np.pi/2))
    # ax.set_xlabel('w[1]', fontsize=30, labelpad=30)
    # ax.set_ylabel('w[2]', fontsize=30, labelpad=30)
    # ax.set_zlabel('w[3]', fontsize=30, labelpad=30)
    ax.set_xlabel('w[value]', fontsize=30, labelpad=30)
    ax.set_ylabel('w[room]', fontsize=30, labelpad=30)
    ax.set_zlabel('w[location]', fontsize=30, labelpad=30)
    ax.tick_params('x', labelsize=25)
    ax.tick_params('y', labelsize=25)
    ax.tick_params('z', labelsize=25)
    ax.set(xlim=(0, 1.03), ylim=(0, 1.03), zlim=(0, 1.03))
    mycm= cm.get_cmap('Reds', 256)
    hmin=min(heats)
    hmax=max(heats)
    print(cells[np.argmin(heats)].bounds)
    print(cells[np.argmax(heats)].bounds)

    cHeats=[(h-hmin)/(hmax-hmin) for h in heats]
    for c, cH in zip(cells, cHeats):
        ax.add_collection3d(Poly3DCollection([c.vertexes], facecolors=mycm(cH), edgecolors='black'))
    # plt.show()
    plt.savefig('./hotel/'+fn + '_k10.' + 'pdf', format='pdf', bbox_inches='tight')


def draw2d_heat(cells, heats, vmin=None, vmax=None, f=None, pdt=None):
    """
    heat='utility' or heat='rkskyband' or heat='rdominate', the default is heat='utility'
    """

    if vmin is None:
        vmin = min(heats)
    if vmax is None:
        vmax = max(heats)
    normH = [(h - vmin) / (vmax - vmin) for h in heats]

    fig, ax = plt.subplots()
    ax.set_xlabel('w[1]')
    ax.set_ylabel('w[2]')
    ax.set(xlim=(0, 1), ylim=(0, 1))

    ax.set_aspect('equal', 'box')
    mycm = cm.get_cmap('Reds', 256)
    for cell, h in zip(cells, normH):
        x = [v[0] for v in cell.vertexes]
        y = [v[1] for v in cell.vertexes]
        ax.plot(x, y, color=mycm(h))
    if pdt is not None:
        ax.plot([i[0] for i in pdt], [i[1] for i in pdt], 'x', markersize=2)
    if not f:
        plt.show()
    else:
        plt.savefig(f+'.png', format='png', bbox_inches='tight')
        # draw2d_heat(cells, heats, vmin=vmin, vmax=vmax)

def draw3d_heat(cells, heats, pdt, fn=None, s=32):
    fig = plt.figure(figsize=(16, 16))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_box_aspect([1, 1, 1])
    ax.view_init(azim=45, elev=90*np.arccos(1.0/(3.0**0.5))/(np.pi/2))
    # ax.set_xlabel('w[1]', fontsize=30, labelpad=30)
    # ax.set_ylabel('w[2]', fontsize=30, labelpad=30)
    # ax.set_zlabel('w[3]', fontsize=30, labelpad=30)
    ax.set_xlabel('w[rebounds]', fontsize=30, labelpad=30)
    ax.set_ylabel('w[assists]', fontsize=30, labelpad=30)
    ax.set_zlabel('w[points]', fontsize=30, labelpad=30)
    ax.tick_params('x', labelsize=25)
    ax.tick_params('y', labelsize=25)
    ax.tick_params('z', labelsize=25)
    ax.set(xlim=(0, 1.03), ylim=(0, 1.03), zlim=(0, 1.03))
    mycm= cm.get_cmap('Reds', 256)
    hmin=min(heats)
    hmax=max(heats)
    # hmin=0.18482252141982863
    # hmax=0.7153214960295415  # TODO remove me, after case study of NBA
    print('heats:', hmin, hmax)
    cHeats=[(h-hmin)/(hmax-hmin) for h in heats]
    for c, cH in zip(cells, cHeats):
        ax.add_collection3d(Poly3DCollection([c.vertexes], facecolors=mycm(cH), edgecolors='black'))
    for p in pdt:
        ax.scatter(p[0], p[1], p[2], color='#069AF3', s=s, marker='x')
    if not fn:
        plt.show()
    else:
        plt.savefig(fn+'.pdf', format='pdf', bbox_inches='tight')

def my2Dheatmap(cells, heats, fn='1'):
    fig = plt.figure(figsize=(16, 16))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_box_aspect([1, 1, 1])
    ax.view_init(azim=45, elev=90*np.arccos(1.0/(3.0**0.5))/(np.pi/2))
    ax.set_xlabel('w[1]', fontsize=30, labelpad=30)
    ax.set_ylabel('w[2]', fontsize=30, labelpad=30)
    ax.set_zlabel('w[3]', fontsize=30, labelpad=30)
    ax.tick_params('x', labelsize=25)
    ax.tick_params('y', labelsize=25)
    ax.tick_params('z', labelsize=25)
    ax.set(xlim=(0, 1.03), ylim=(0, 1.03), zlim=(0, 1.03))
    mycm= cm.get_cmap('Reds', 256)
    hmin=min(heats)
    hmax=max(heats)
    cHeats=[(h-hmin)/(hmax-hmin) for h in heats]
    for c, cH in zip(cells, cHeats):
        ax.add_collection3d(Poly3DCollection([c.vertexes], facecolors=mycm(cH), edgecolors='black'))
    # plt.show()
    plt.savefig('./cs1_3d/'+fn + '_k10.' + 'png', format='png', bbox_inches='tight')


def my3Dheatmap2(cells, heats, pdts, colors, markers, fn='1', s=32):
    fig = plt.figure(figsize=(16, 16))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_box_aspect([1, 1, 1])
    ax.view_init(azim=45, elev=90*np.arccos(1.0/(3.0**0.5))/(np.pi/2))
    ax.set_xlabel('w[1]', fontsize=30, labelpad=30)
    ax.set_ylabel('w[2]', fontsize=30, labelpad=30)
    ax.set_zlabel('w[3]', fontsize=30, labelpad=30)
    ax.tick_params('x', labelsize=25)
    ax.tick_params('y', labelsize=25)
    ax.tick_params('z', labelsize=25)
    ax.set(xlim=(0, 1.03), ylim=(0, 1.03), zlim=(0, 1.03))
    mycm= cm.get_cmap('Reds', 256)
    hmin=min(heats)
    hmax=max(heats)
    cHeats=[(h-hmin)/(hmax-hmin) for h in heats]
    for c, cH in zip(cells, cHeats):
        ax.add_collection3d(Poly3DCollection([c.vertexes], facecolors=mycm(cH), edgecolors='black'))
    for pdt, color, marker in zip(pdts, colors, markers):
        pdt_t=np.array(pdt).T
        ax.scatter(pdt_t[0], pdt_t[1], pdt_t[2], color=color, s=s, marker=marker)
    plt.show()
    # plt.savefig('./cs2_3d/'+fn + '.' + 'pdf', format='pdf', bbox_inches='tight')


def inCell(cell, u):
    s=sum(u)
    u=[i/s for i in u]
    for d in range(cell.dim):
        if not (cell.bounds[d*2] <= u[d] < cell.bounds[d*2+1]):
            return False
    return True



def count_user(cells, user):
    ret=[0 for _ in range(len(cells))]
    for u in user:
        # if sum(u)<0.5:
        #     continue
        for c, i in zip(cells, range(len(ret))):
            if inCell(c, u):
                ret[i]+=1
                break
    return ret

def inCell2(parent, u):
    for child in parent.children:
        if inCell(child, u):
            if child.cur_l>=child.tar_l:
                child.cnt+=1
            else:
                inCell2(child, u)
            break


def count_user2(root, user):
    for u in user:
        inCell2(root, u)
import random
def draw_user(user, sample=0.01, f=None):
    norm_user=[]
    for u in user:
        s=sum(u)
        norm_user.append([i/s for i in u])
    new_usr=[]
    if sample:
        random.seed(1)
        random.shuffle(norm_user)
        for i in range(int(len(norm_user)*sample)):
            new_usr.append(norm_user[i])
    else:
        new_usr=norm_user
    if f:
        draw_3d_pdt(new_usr, f=f+'.pdf')
        draw_3d_pdt(new_usr, f=f+'_2.pdf', view=True)
    else:
        draw_3d_pdt(new_usr, f=f)
        draw_3d_pdt(new_usr, f=f, view=True)

def cs1_3dplot(cells, pdt, user, attrs='0_1_2', root=None):
    uh=[c.uHeat for c in cells]
    rh=[c.rHeat for c in cells]
    ch=[c.dHeat for c in cells]
    # print('begin user count')
    # userH=count_user(cells, user)
    norm_user=[]
    for u in user:
        s=sum(u)
        norm_user.append([i/s for i in u])
    count_user2(root, norm_user)
    userH=[c.cnt for c in cells]
    statist=[0 for i in range(len(norm_user[0]))]
    for u in norm_user:
        idx=0
        maxv=0
        for i in range(len(u)):
            if maxv<u[i]:
                idx=i
                maxv=u[i]
        statist[idx]+=1
    print('statist:', statist)


    print('$', np.corrcoef([uh, rh]))
    print('$', np.corrcoef([uh, ch]))
    print('$', np.corrcoef([uh, userH]))
    print('$', np.corrcoef([rh, ch]))
    print('$', np.corrcoef([rh, userH]))
    print('$', np.corrcoef([ch, userH]))
    # distance.jensenshannon is to calculate JS distance, so there should be a "**2" so as to get JS divergence
    # to normalize the JS divergence from 0 to 1, we choose base 2
    print(distance.jensenshannon([1.0 for _ in range(len(userH))], userH, base=2)**2, bhattacharyya_coefficient([1.0 for _ in range(len(userH))], userH))

    min_uh=min(uh)
    max_uh=max(uh)
    uh=[(i-min_uh)/(max_uh-min_uh) for i in uh]
    min_ch=min(ch)
    max_ch=max(ch)
    ch=[(i-min_ch)/(max_ch-min_ch) for i in ch]
    min_rh=min(rh)
    max_rh=max(rh)
    rh=[(i-min_rh)/(max_rh-min_rh) for i in rh]
    print([min(uh), max(uh)], distance.jensenshannon(uh, userH, base=2)**2, bhattacharyya_coefficient(uh, userH))
    print([min(rh), max(rh)], distance.jensenshannon(rh, userH, base=2)**2, bhattacharyya_coefficient(rh, userH))
    print([min(ch), max(ch)], distance.jensenshannon(ch, userH, base=2)**2, bhattacharyya_coefficient(ch, userH))
    # uh_min_pdt=pdt[cells[np.argmin(uh)].rkskyband]
    # uh_max_pdt=pdt[cells[np.argmax(uh)].rkskyband]
    # rh_min_pdt=pdt[cells[np.argmin(rh)].rkskyband]
    # rh_max_pdt=pdt[cells[np.argmax(rh)].rkskyband]
    # ch_min_pdt=pdt[cells[np.argmin(ch)].rkskyband]
    # ch_max_pdt=pdt[cells[np.argmax(ch)].rkskyband]
    # print(np.corrcoef(uh_min_pdt.T), np.corrcoef(uh_max_pdt.T), sep='\r\n', end='*'*50+'\n')
    # print(np.corrcoef(rh_min_pdt.T), np.corrcoef(rh_max_pdt.T), sep='\r\n', end='*'*50+'\n')
    # print(np.corrcoef(ch_min_pdt.T), np.corrcoef(ch_max_pdt.T), sep='\r\n', end='*'*50+'\n')



    # print([min(userH), max(userH)], distance.jensenshannon(userH, userH))

    # draw_3d_pdt(norm_user, f='./hotel/user_n.pdf')
    # draw_3d_pdt(norm_user, f='./hotel/user2_n.pdf', view=True)
    #
    # # TODO delete me later
    my3Dheatmap(cells, uh, attrs+'uh')
    #
    my3Dheatmap(cells, rh, attrs+'rh')
    my3Dheatmap(cells, ch, attrs+'ch')
    # my3Dheatmap(cells, ch, attrs+'ch')

    # my3Dheatmap(cells, userH, attrs+'userH')


def cs1_plot(cells, pdt, user, attrs='0_1_2', root=None):
    if len(pdt)==0:
        return
    if len(pdt[0])==2: # 2d-plot
        pass
    elif len(pdt[0]==3): # 3d-plot
        cs1_3dplot(cells, pdt, user, attrs, root=root)

def cs2_3dplot(cells, pdt, user, attrs='0_1_2'):
    uh=[c.uHeat for c in cells]
    rh=[c.rHeat for c in cells]
    ch=[c.dHeat for c in cells]
    userH=count_user(cells, user)
    uh_min_pdt=pdt[cells[np.argmin(uh)].rkskyband]
    uh_max_pdt=pdt[cells[np.argmax(uh)].rkskyband]
    rh_min_pdt=pdt[cells[np.argmin(rh)].rkskyband]
    rh_max_pdt=pdt[cells[np.argmax(rh)].rkskyband]
    ch_min_pdt=pdt[cells[np.argmin(ch)].rkskyband]
    ch_max_pdt=pdt[cells[np.argmax(ch)].rkskyband]
    print(min(userH), max(userH))
    print(np.corrcoef(uh_min_pdt.T), np.corrcoef(uh_max_pdt.T), sep='\r\n')
    print(np.corrcoef(rh_min_pdt.T), np.corrcoef(rh_max_pdt.T), sep='\r\n')
    print(np.corrcoef(ch_min_pdt.T), np.corrcoef(ch_max_pdt.T), sep='\r\n')
    print('$', np.corrcoef([rh, ch]))

    my3Dheatmap2(cells, uh, [uh_min_pdt, uh_max_pdt], ['blue', 'red'], [4, 5],  s=128, fn=attrs+'uh')
    my3Dheatmap2(cells, rh, [rh_min_pdt, rh_max_pdt], ['blue', 'red'], [4, 5],  s=128, fn=attrs+'rh')
    my3Dheatmap2(cells, ch, [ch_min_pdt, ch_max_pdt], ['blue', 'red'], [4, 5],  s=128, fn=attrs+'ch')


def cs2_plot(cells, pdt, user, attrs='0_1_2'):
    if len(pdt)==0:
        return
    if len(pdt[0])==2: # 2d-plot
        pass
    elif len(pdt[0]==3): # 3d-plot
        cs2_3dplot(cells, pdt, user, attrs)

def cs_2d_5cluster(cells, pdt=None, f=None, norm=True):
    uh=[c.uHeat for c in cells]
    rh=[c.rHeat for c in cells]
    ch=[c.dHeat for c in cells]
    if not f:
        draw2d_heat(cells, uh, pdt=pdt)
        draw2d_heat(cells, rh, pdt=pdt)
        draw2d_heat(cells, ch, pdt=pdt)
    else:
        if norm:
            draw2d_heat(cells, uh, pdt=pdt, f=f+"_uh_n")
        else:
            draw2d_heat(cells, uh, pdt=pdt, f=f+"_uh")
        draw2d_heat(cells, rh, pdt=pdt, f=f+"_rh")
        draw2d_heat(cells, ch, pdt=pdt, f=f+"_ch")


def cs_3d_5cluster(cells, pdt=None, f=None, norm=True):
    uh=[c.uHeat for c in cells]
    rh=[c.rHeat for c in cells]
    ch=[c.dHeat for c in cells]
    if not f:
        draw3d_heat(cells, uh, pdt=pdt)
        draw3d_heat(cells, rh, pdt=pdt)
        draw3d_heat(cells, ch, pdt=pdt)
    else:
        if norm:
            draw3d_heat(cells, uh, pdt=pdt, fn=f+"_uh_n")
        else:
            draw3d_heat(cells, uh, pdt=pdt, fn=f+"_uh")
        draw3d_heat(cells, rh, pdt=pdt, fn=f+"_rh")
        draw3d_heat(cells, ch, pdt=pdt, fn=f+"_ch")

def draw_2d_pdt(pdt_ta, f=None):
    fig = plt.figure()
    ax = fig.gca()
    ax.plot([i[0] for i in pdt_ta], [i[1] for i in pdt_ta], 'x', markersize=2)
    ax.axis((0, 1, 0, 1))
    ax.set_aspect('equal')
    if not f:
        plt.show()
    else:
        plt.savefig(f, format='png', bbox_inches='tight')

def draw_3d_pdt(pdt_ta, f=None, view=False):
    fig = plt.figure(figsize=(16, 16))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_box_aspect([1, 1, 1])
    if view:
        ax.view_init(azim=45, elev=90 * np.arccos(1.0 / (3.0 ** 0.5)) / (np.pi / 2))
    # ax.set_xlabel('w[1]', fontsize=30, labelpad=30)
    # ax.set_ylabel('w[2]', fontsize=30, labelpad=30)
    # ax.set_zlabel('w[3]', fontsize=30, labelpad=30)
    ax.set_xlabel('w[value]', fontsize=30, labelpad=30)
    ax.set_ylabel('w[room]', fontsize=30, labelpad=30)
    ax.set_zlabel('w[location]', fontsize=30, labelpad=30)
    ax.tick_params('x', labelsize=25)
    ax.tick_params('y', labelsize=25)
    ax.tick_params('z', labelsize=25)
    ax.set(xlim=(0, 1.03), ylim=(0, 1.03), zlim=(0, 1.03))
    ax.plot([i[0] for i in pdt_ta], [i[1] for i in pdt_ta], [i[2] for i in pdt_ta], 'x')
    if not f:
        plt.show()
    else:
        plt.savefig(f, format='pdf', bbox_inches='tight')
