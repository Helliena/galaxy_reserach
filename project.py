import numpy as np
import json
import astropy
from timeit import default_timer as t
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import WMAP5, WMAP7
import matplotlib.pyplot as plt
import copy 
import random


# n = int(input('Введите количество кораблей: '))
n = 100
speed_of_light = 300_000
fraction_speed_of_light = 0.01
year = 365*24*3600


# мы берем звезды и фигачим их данные в 4 массива
ra = []
s = t()
dec = []
source_id = []
par = []

l = []
b = []


with open('gaiadr3_gaia_source.json') as f:
    templates = json.load(f)
v0 = templates.get("data")

print(len(v0))


'''
Нужно создать данные о 5000-10000 ближайших звёзд

'''


# Добавление только звёзд с полными данными
for i in range(0, 1000, 1):
    if v0[i][3] != None:
        source_id.append(v0[i][0])
        ra.append(v0[i][1])
        dec.append(v0[i][2])
        par.append(v0[i][3])


# Перевод от экваториальных к галактическим координатам
for i in range(len(dec)):
    c_icrs = SkyCoord(ra=ra[i]*u.degree, dec=dec[i]*u.degree, frame='icrs')
    l.append(c_icrs.galactic.l.degree)
    b.append(c_icrs.galactic.b.degree)
    
l = np.array(l)
b = np.array(b) 
par = np.array(par)/1000
dist = 1/par * 3.26 # расстояние от солнца до звёзд в световых годах

# Проецируем на плоскость и переходим от полярных координат к декартовым
dist_proj = np.cos(b) * dist
x = dist_proj * np.cos(l)
y = dist_proj * np.sin(l)



minim = np.zeros((7,len(dec)))



c1 = SkyCoord(ra=0*u.degree, dec=0*u.degree, distance=0*u.pc, frame='icrs')
c2 = SkyCoord(ra=292.674*u.degree, dec=dec[i]*u.degree, distance=dist[i]*u.pc, frame='icrs')
r = c1.separation_3d(c2)
f = r.pc

# Рассчитали расстояния от звёзд до Солнца
for i in range(len(dec)):
    c1 = SkyCoord(ra=0*u.degree, dec=0*u.degree, distance=0*u.pc, frame='icrs')
    c2 = SkyCoord(ra=ra[i]*u.degree, dec=dec[i]*u.degree, distance=dist[i]*u.pc, frame='icrs')
    r = c1.separation_3d(c2)
    f = r.pc
    minim[0][i] = f
    minim[1][i] = source_id[i]
    minim[2][i] = ra[i]
    minim[3][i] = dec[i]
    minim[4][i] = dist[i]
    minim[5][i] = x[i]
    minim[6][i] = y[i]

# minim1 - список со всеми звёздами из которого потом подгружаются данные
# Сортировка по возрастанию расстояний от солнца до звёзд с запоминанием id каждой звезды     
minim1 = minim[:, minim[0].argsort()]

# Для проверки id звезды требуется использовать int
# print(int(minim[1][i]) == 5959655667469483003)
    
# Вероятность
p = 0.3

id_defeat = []
id_success = []

# время передвижения корабля 
travel_time = []


# добавляем в списки id звёзд успешных и неудачных
for i in range(n):  
    a = random.random()
    if a > p:
        id_defeat.append(minim1[1][i])
    else:
        travel_time.append(minim1[0][i]) # время в пути в годах
        id_success.append(minim1[1][i])

    
# Удаляем звёзды которые неудачно колонизировали
for i in id_defeat:
    if int(i) in minim[1]:  
        minim = np.delete(minim, [np.where(minim == int(i))[1]], axis = 1)

    
# Удаляем звёзды которые колонизировали
for i in id_success:
    if int(i) in minim[1]:  
        minim = np.delete(minim, [np.where(minim == int(i))[1]], axis = 1) 
    
# Уменьшеаем число кораблей
n = n - len(id_defeat) 


id_success_new = []
id_defeat_new = []


id_success_all = id_success.copy()
k = 0
while n>0:
    print(n)
    for i in id_success:
        # берём значения успешных звёзд из целого массива 
        index = np.where(minim1 == int(i))[1]
        new_star = np.zeros((2,len(minim[1])))

        for j in range(len(minim[1])):
            # координаты успешных звёзд
            c3 = SkyCoord(ra=minim1[2][index]*u.degree, dec=minim1[3][index]*u.degree, distance=minim1[4][index]*u.pc, frame='icrs')
            # координаты оставшихся звёзд
            c4 = SkyCoord(ra=minim[2][j]*u.degree, dec=minim[3][j]*u.degree, distance=minim[4][j]*u.pc, frame='icrs')
            
            o = c3.separation_3d(c4)
            v = o.pc # рассчитатываем новые ближние звёзды для каждой успешной звезды

            # Записываем новые ближайшие звёзды и сортируем
            new_star[0][j] = v[0]
            new_star[1][j] = minim[1][j]
            new_star1 = new_star[:, new_star[0].argsort()]
            
            
        # неудачные звёзды, добавляем в отдельный список удачные
        a = random.random()
        # print(a)
        if a > p:
 
            id_defeat_new.append(new_star1[1][0])
        else:
            travel_time.append(travel_time[k]+new_star1[0][0])
            # print(new_star1[0][0])
            id_success_new.append(new_star1[1][0])
            id_success_all.append(new_star1[1][0])
        k+=1 
        
        # print(id_success_new)
        for i in id_defeat_new:
            if int(i) in minim[1]:  
                minim = np.delete(minim, [np.where(minim == int(i))[1]], axis = 1)
        for i in id_success_new:
            if int(i) in minim[1]:  
                minim = np.delete(minim, [np.where(minim == int(i))[1]], axis = 1)
                
                
    id_success = id_success_new
    id_success_new = []
    n = n - len(id_defeat_new) 
    id_defeat_new = []
    
    
    
print(n)




'''
Для отображения требуется соединить в один массив id всех успешных кораблей, 
их x и y координаты и travel_time. Сортировать всё по travel_time и построить гифку
'''
plt.style.use(['dark_background'])
fig = plt.figure(figsize=(50, 50), frameon = True)
plt.scatter(x, y, color='white', s=100)
plt.scatter(0, 0, color='yellow', s=125)
plt.xticks(fontsize = 50, rotation=30)
plt.yticks(fontsize = 50)

# plt.scatter(8827 * 3.26, 0, color='blue', s=5)
# plt.plot(x[-1], y[-1], color = 'blue', marker = 'o')
plt.xlim(-10*3.26, 10*3.26)
plt.ylim(-10*3.26, 10*3.26)
plt.text(-6,8*3.26, f'0 лет', size=100)
plt.title(f'Вероятность = {p}', size=100)
plt.savefig(f'0.png')

fig.clear()
plt.close(fig)


data = np.zeros((4,len(id_success_all)))

for i in range(len(id_success_all)):
    index = np.where(minim1 == int(id_success_all[i]))[1]
    data[1,i] = id_success_all[i]  
    data[0,i] = travel_time[i]  
    data[2,i] = minim1[5][index][0]
    data[3,i] = minim1[6][index][0]

data = data[:, data[0].argsort()]

for i in range(len(id_success_all)):
    index = np.where(minim1 == int(id_success_all[i]))[1]
    # data[2,i] = minim1[5][index][0]
    # data[3,i] = minim1[6][index][0]
    fig = plt.figure(figsize=(50,50), frameon = True)
    plt.scatter(x, y, color='white', s=100)
    plt.scatter(0, 0, color='yellow', s=125)
    plt.xticks(fontsize = 50, rotation=30)
    plt.yticks(fontsize = 50)
    # plt.scatter(8827 * 3.215006, 0, color='blue', s=5)
    plt.scatter(data[2], data[3], color='red', s=100)
    plt.text(-6,8*3.26, f'{round(data[0,i]/fraction_speed_of_light)} лет', size=100)
    plt.xlim(-10*3.26, 10*3.26)
    plt.ylim(-10*3.26, 10*3.26)
    plt.title(f'Вероятность = {p}', size=100)
    plt.savefig(f'{i+1}.png')
    fig.clear()
    plt.close(fig)
    plt.pause(0.1)
  
    


d = t() - s
print(f'Время работы программы: {d} секунд')