import math
import numpy as np
import matplotlib.pyplot as plt
#常量
acc_scale=1.5258789063e-06
gyo_scale=1.0850694444e-07
f=100
g=9.7936174
we=math.degrees(7.292115e-5)
fai=math.radians(30.531651244)
wes=we*math.sin(fai)

#静态测量的加速度计数据
class static_f:
    def __init__(self):
        self.fx = 0
        self.fy = 0
        self.fz = 0

#静态测量的陀螺数据
class static_g:
    def __init__(self):
        self.wx = 0
        self.wy = 0
        self.wz = 0

#动态测量的陀螺数据
class dynamic_g:
    def __init__(self):
        self.wx= 0
        self.wy= 0
        self.wz= 0

#读取静态数据文件
def read_staticdata(filepath, static_f_list,static_g_list):
    with open(filepath, 'r') as file:
        for line in file:
            if line.startswith("%RAWIMUSA"):

                obf=static_f()
                obg=static_g()

                line_data = line.strip().split(',')
                #读入为比力，单位m/s2
                obf.fx = float(line_data[7]) * acc_scale * f
                obf.fy = -float(line_data[6]) * acc_scale * f
                obf.fz = float(line_data[5]) * acc_scale * f
                static_f_list.append(obf)
                #读入角速度，单位deg/s
                stri = line_data[10].split("*")
                obg.wx=math.degrees(float(stri[0]) * gyo_scale * f)
                obg.wy=math.degrees(-float(line_data[9]) * gyo_scale * f)
                obg.wz=math.degrees(float(line_data[8]) * gyo_scale * f)
                static_g_list.append(obg)

#读取动态（转动）数据文件
def read_dynamicdata(filepath, dynamic_g_list):
    with open(filepath, 'r') as file:
        for line in file:
            if line.startswith("%RAWIMUSA"):
                obg=dynamic_g()
                line_data = line.strip().split(',')
                stri = line_data[10].split("*")
                # 读入角速度，单位deg/s
                obg.wx=math.degrees(float(stri[0]) * gyo_scale * f)
                obg.wy=math.degrees(-float(line_data[9]) * gyo_scale * f)
                obg.wz=math.degrees(float(line_data[8]) * gyo_scale * f)
                dynamic_g_list.append(obg)

#求列表元素均值
def mean(list):
    m = sum(list) / len(list)
    return m

#六位置法标定加速度计
def standard_acc(xup_list,xdown_list,yup_list,ydown_list,zup_list,zdown_list):
    xup_fx = [xup.fx for xup in xup_list]
    #print(sum(xup_fx),len(xup_fx))
    xup_fy = [xup.fy for xup in xup_list]
    xup_fz = [xup.fz for xup in xup_list]
    xdown_fx = [xdown.fx for xdown in xdown_list]
    xdown_fy = [xdown.fy for xdown in xdown_list]
    xdown_fz = [xdown.fz for xdown in xdown_list]
    yup_fx = [yup.fx for yup in yup_list]
    yup_fy = [yup.fy for yup in yup_list]
    yup_fz = [yup.fz for yup in yup_list]
    ydown_fx = [ydown.fx for ydown in ydown_list]
    ydown_fy = [ydown.fy for ydown in ydown_list]
    ydown_fz = [ydown.fz for ydown in ydown_list]
    zup_fx = [zup.fx for zup in zup_list]
    zup_fy = [zup.fy for zup in zup_list]
    zup_fz = [zup.fz for zup in zup_list]
    zdown_fx = [zdown.fx for zdown in zdown_list]
    zdown_fy = [zdown.fy for zdown in zdown_list]
    zdown_fz = [zdown.fz for zdown in zdown_list]

    A = np.array([[g, -g,0,0,0,0], [0,0,g,-g,0,0],[0,0,0,0,g,-g],[1,1,1,1,1,1]])
    L = np.array([[mean(xup_fx), mean(xdown_fx), mean(yup_fx), mean(ydown_fx), mean(zup_fx), mean(zdown_fx)],
                  [mean(xup_fy), mean(xdown_fy), mean(yup_fy), mean(ydown_fy), mean(zup_fy), mean(zdown_fy)],
                  [mean(xup_fz), mean(xdown_fz), mean(yup_fz), mean(ydown_fz), mean(zup_fz), mean(zdown_fz)]])

    LAT=np.dot(L, A.transpose())
    AAT=np.dot(A, A.transpose())
    M=np.dot(LAT, np.linalg.inv(AAT))
    #按该顺序输出：bx,by,bz,sx,sy,sz,ryx,rzx,rxy,rzy,rxz,ryz
    #print(M)
    return M[0,3],M[1,3],M[2,3],M[0,0]-1,M[1,1]-1,M[2,2]-1,M[0,1],M[0,2],M[1,0],M[1,2],M[2,0],M[2,1]

#两位置法标定陀螺零偏
def standard_gyroB(xup_list,xdown_list,yup_list,ydown_list,zup_list,zdown_list):
    xup_wx = [xup.wx for xup in xup_list]
    xdown_wx = [xdown.wx for xdown in xdown_list]
    yup_wy = [yup.wy for yup in yup_list]
    ydown_wy = [ydown.wy for ydown in ydown_list]
    zup_wz = [zup.wz for zup in zup_list]
    zdown_wz = [zdown.wz for zdown in zdown_list]
    bx = (mean(xup_wx) + mean(xdown_wx)) / 2
    by = (mean(yup_wy) + mean(ydown_wy)) / 2
    bz = (mean(zup_wz) + mean(zdown_wz)) / 2
    #print(math.radians(bx),math.radians(by),math.radians(bz))
    return bx,by,bz

#获取转动开始和结束的索引下标
def ifrevolve(list):
    i=0
    start=0
    end=0
    for value in list:
        if(start==0):
            if(abs(value)>=0.00001):
                start=i
        if(end==0):
            if(abs(value)<=0.00001 and (i-start)>3000):
                end=i
        i=i+1
    return start,end

#角位置法辅助函数
def help_standard_gyro(cw_list,ccw_list):
    cw_s, cw_e = ifrevolve(cw_list)
    ccw_s, ccw_e = ifrevolve(ccw_list)
    #print(cw_s, cw_e,ccw_s, ccw_e)
    # 保证数据量相同
    if ((cw_e - cw_s) < (ccw_e - ccw_s)):
        cw_e = cw_e + (ccw_e - ccw_s) - (cw_e - cw_s)
    if ((cw_e - cw_s) > (ccw_e - ccw_s)):
        ccw_e = ccw_e + (cw_e - cw_s) - (ccw_e - ccw_s)
    valid_cw = cw_list[cw_s:cw_e + 1]
    valid_ccw = ccw_list[ccw_s: ccw_e + 1]
    t=(len(valid_cw)-1)/f
    #t1=(len(valid_ccw)-1)/f
    b = (sum(valid_cw) + sum(valid_ccw)) / (2*f*t) - wes
    s = (sum(valid_cw) - sum(valid_ccw)) / (2 * f * 360)-1
    #print(sum(valid_cw)/f,sum(valid_ccw)/f)
    #print(t,t1,wes)
    #print((sum(valid_cw)/f,"+",sum(valid_ccw)/f,"/",2*t,"-",wes))
    return b,s

#角位置法标定陀螺零偏、比例因子误差
def standard_gyro(xcw_list,xccw_list,ycw_list,yccw_list,zcw_list,zccw_list):
    xcw_wx = [xcw.wx for xcw in xcw_list]
    xccw_wx = [xccw.wx for xccw in xccw_list]
    ycw_wy = [ycw.wy for ycw in ycw_list]
    yccw_wy = [yccw.wy for yccw in yccw_list]
    zcw_wz = [zcw.wz for zcw in zcw_list]
    zccw_wz = [zccw.wz for zccw in zccw_list]
    bx, sx = help_standard_gyro(xcw_wx,xccw_wx)
    by, sy = help_standard_gyro(ycw_wy, yccw_wy)
    bz, sz = help_standard_gyro(zcw_wz,zccw_wz)
    return bx,by,bz,sx,sy,sz

#加速度计的误差补偿
def compensate_acc(original_acc,new_acc,bx,by,bz,sx,sy,sz,ryx,rzx,rxy,rzy,rxz,ryz):
    M = np.array([[1+sx,ryx,rzx],
                 [rxy, 1+sy, rzy],
                 [rxz, ryz, 1+sz]])
    M_=np.linalg.inv(M)
    B = np.array([[bx],[by],[bz]])
    #print(B)
    for i in range(len(original_acc)):
        ori=np.array([[original_acc[i].fx],[original_acc[i].fy],[original_acc[i].fz]])
        new=np.dot(M_, ori-B)
        obj=static_f()
        obj.fx=new[0,0]
        obj.fy=new[1,0]
        obj.fz=new[2,0]
        #print(obj)
        new_acc.append(obj)

#陀螺的误差补偿
def compensate_gyro(original_gyro,new_gyro,bx,by,bz,sx,sy,sz):
    for i in range(len(original_gyro)):
        obj=dynamic_g()
        obj.wx=(original_gyro[i].wx-bx)/(1+sx)
        obj.wy=(original_gyro[i].wy-by)/(1+sy)
        obj.wz=(original_gyro[i].wz-bz)/(1+sz)
        new_gyro.append(obj)


filepath1 = r'D:\惯性导航编程作业\实验一 标定\data\x_up_static.ASC'
acc_x_up=[]
static_g_x_up=[]
read_staticdata(filepath1, acc_x_up,static_g_x_up)

filepath2 = r'D:\惯性导航编程作业\实验一 标定\data\x_down_static.ASC'
acc_x_down=[]
static_g_x_down=[]
read_staticdata(filepath2, acc_x_down,static_g_x_down)

filepath3 = r'D:\惯性导航编程作业\实验一 标定\data\y_up_static.ASC'
acc_y_up=[]
static_g_y_up=[]
read_staticdata(filepath3, acc_y_up,static_g_y_up)

filepath4 = r'D:\惯性导航编程作业\实验一 标定\data\y_down_static.ASC'
acc_y_down=[]
static_g_y_down=[]
read_staticdata(filepath4, acc_y_down,static_g_y_down)

filepath5 = r'D:\惯性导航编程作业\实验一 标定\data\z_up_static.ASC'
acc_z_up=[]
static_g_z_up=[]
read_staticdata(filepath5, acc_z_up,static_g_z_up)

filepath6 = r'D:\惯性导航编程作业\实验一 标定\data\z_down_static.ASC'
acc_z_down=[]
static_g_z_down=[]
read_staticdata(filepath6, acc_z_down,static_g_z_down)

bax,bay,baz,sax,say,saz,rayx,razx,raxy,razy,raxz,rayz=standard_acc(acc_x_up,acc_x_down,acc_y_up,acc_y_down,acc_z_up,acc_z_down)
print("加速度计的bax(m/s2),bay(m/s2),baz(m/s2),sax,say,saz,rayx,razx,raxy,razy,raxz,rayz分别为：")
print("{:.9f} {:.9f} {:.9f} {:.9f} {:.9f} {:.9f} {:.9f} {:.9f} {:.9f} {:.9f} {:.9f} {:.9f}".format(bax, bay, baz, sax, say, saz, rayx, razx, raxy, razy, raxz, rayz))

bgx1,bgy1,bgz1=standard_gyroB(static_g_x_up,static_g_x_down,static_g_y_up,static_g_y_down,static_g_z_up,static_g_z_down)
print("两位置法标定的陀螺的三轴零偏分别为(bgx(deg/h),bgy(deg/h),bgz(deg/h))：")
print("{:.9f} {:.9f} {:.9f}".format(bgx1*3600, bgy1*3600, bgz1*3600))


filepath7 = r'D:\惯性导航编程作业\实验一 标定\data\x_+.ASC'
dynamic_g_xcw=[]
read_dynamicdata(filepath7, dynamic_g_xcw)

filepath8 = r'D:\惯性导航编程作业\实验一 标定\data\x_-.ASC'
dynamic_g_xccw=[]
read_dynamicdata(filepath8, dynamic_g_xccw)

filepath9 = r'D:\惯性导航编程作业\实验一 标定\data\y_+.ASC'
dynamic_g_ycw=[]
read_dynamicdata(filepath9, dynamic_g_ycw)

filepath10 = r'D:\惯性导航编程作业\实验一 标定\data\y_-.ASC'
dynamic_g_yccw=[]
read_dynamicdata(filepath10, dynamic_g_yccw)

filepath11 = r'D:\惯性导航编程作业\实验一 标定\data\z_+.ASC'
dynamic_g_zcw=[]
read_dynamicdata(filepath11, dynamic_g_zcw)

filepath12 = r'D:\惯性导航编程作业\实验一 标定\data\z_-.ASC'
dynamic_g_zccw=[]
read_dynamicdata(filepath12, dynamic_g_zccw)

bgx,bgy,bgz,sgx,sgy,sgz=standard_gyro(dynamic_g_xcw,dynamic_g_xccw,dynamic_g_ycw,dynamic_g_yccw,dynamic_g_zcw,dynamic_g_zccw)
print("角位置法标定的陀螺三轴零偏和比例因子分别为（bgx(deg/h),bgy(deg/h),bgz(deg/h),sgx,sgy,sgz）：")
print("{:.9f} {:.9f} {:.9f} {:.9f} {:.9f} {:.9f}".format(bgx*3600, bgy*3600, bgz*3600, sgx, sgy, sgz))

'''
#画出z轴朝上的正转陀螺输出
time=[]
for t in range(len(dynamic_g_zcw)):
    time.append(t/f)
owx=[value.wx for value in dynamic_g_zcw]
owy=[value.wy for value in dynamic_g_zcw]
owz=[value.wz for value in dynamic_g_zcw]
figure=plt.figure(figsize=(10,5),num="w")  # 设置画布的大小
plt.plot(time, owx, label="wx")
plt.plot(time, owy, label="wy")
plt.plot(time, owz, label="wz")
# 添加标题和标签
plt.title("w")  # 设置整个图的标题
plt.xlabel("t(s)")  # x 轴标签
plt.ylabel("w(m/s2)")  # y 轴标签
# 添加图例
plt.legend()
# 显示图形
plt.show()
'''

#补偿加速度计误差并作前后对比图,以x轴朝上的静态数据为例
#acc_x_up=[]
new_acc_x_up=[]
compensate_acc(acc_x_up,new_acc_x_up,bax,bay,baz,sax,say,saz,rayx,razx,raxy,razy,raxz,rayz)
time1=[]
for t in range(len(new_acc_x_up)):
    time1.append(t/f)
ori_fx=[value.fx for value in acc_x_up]
new_fx=[value.fx for value in new_acc_x_up]
ori_fy=[value.fy for value in acc_x_up]
new_fy=[value.fy for value in new_acc_x_up]
ori_fz=[value.fz for value in acc_x_up]
new_fz=[value.fz for value in new_acc_x_up]
# 绘制曲线
#fx
figure1=plt.figure(figsize=(10,5),num="original and newdata after compensate(fx)")  # 设置画布的大小
plt.plot(time1, ori_fx, label="ori_fx")
plt.plot(time1, new_fx, label="new_fx")
# 添加标题和标签
plt.title("original and newdata after compensate(fx)")  # 设置整个图的标题
plt.xlabel("t(s)")  # x 轴标签
plt.ylabel("fx(m/s2)")  # y 轴标签
# 添加图例
plt.legend()
# 显示图形
plt.show()
#fy
figure2=plt.figure(figsize=(10,5),num="original and newdata after compensate(fy)")  # 设置画布的大小
plt.plot(time1, ori_fy, label="ori_fy")
plt.plot(time1, new_fy, label="new_fy")
# 添加标题和标签
plt.title("original and newdata after compensate(fy)")  # 设置整个图的标题
plt.xlabel("t(s)")  # x 轴标签
plt.ylabel("fy(m/s2)")  # y 轴标签
# 添加图例
plt.legend()
# 显示图形
plt.show()
#fz
figure3=plt.figure(figsize=(10,5),num="original and newdata after compensate(fz)")  # 设置画布的大小
plt.plot(time1, ori_fz, label="ori_fz")
plt.plot(time1, new_fz, label="new_fz")
# 添加标题和标签
plt.title("original and newdata after compensate(fz)")  # 设置整个图的标题
plt.xlabel("t(s)")  # x 轴标签
plt.ylabel("fz(m/s2)")  # y 轴标签
# 添加图例
plt.legend()
# 显示图形
plt.show()

#补偿陀螺误差并作前后对比图,以x轴+转的数据为例
#dynamic_g_xcw=[]
new_dynamic_g_xcw=[]
compensate_gyro(dynamic_g_xcw,new_dynamic_g_xcw,bgx1,bgy1,bgz1,sgx,sgy,sgz)
time2=[]
for t in range(len(new_dynamic_g_xcw)):
    time2.append(t/f)
ori_wx=[value.wx for value in dynamic_g_xcw]
new_wx=[value.wx for value in new_dynamic_g_xcw]
ori_wy=[value.wy for value in dynamic_g_xcw]
new_wy=[value.wy for value in new_dynamic_g_xcw]
ori_wz=[value.wz for value in dynamic_g_xcw]
new_wz=[value.wz for value in new_dynamic_g_xcw]
# 绘制曲线
#wx
figure4=plt.figure(figsize=(10,5),num="original and newdata after compensate(wx)")  # 设置画布的大小
plt.plot(time2, ori_wx, label="ori_wx")
plt.plot(time2, new_wx, label="new_wx")
# 添加标题和标签
plt.title("original and newdata after compensate(wx)")  # 设置整个图的标题
plt.xlabel("t(s)")  # x 轴标签
plt.ylabel("wx(deg/s)")  # y 轴标签
# 添加图例
plt.legend()
# 显示图形
plt.show()
#wy
figure5=plt.figure(figsize=(10,5),num="original and newdata after compensate(wy)")  # 设置画布的大小
plt.plot(time2, ori_wy, label="ori_wy")
plt.plot(time2, new_wy, label="new_wy")
# 添加标题和标签
plt.title("original and newdata after compensate(wy)")  # 设置整个图的标题
plt.xlabel("t(s)")  # x 轴标签
plt.ylabel("wy(deg/s)")  # y 轴标签
# 添加图例
plt.legend()
# 显示图形
plt.show()
#wz
figure6=plt.figure(figsize=(10,5),num="original and newdata after compensate(wz)")  # 设置画布的大小
plt.plot(time2, ori_wz, label="ori_wz")
plt.plot(time2, new_wz, label="new_wz")
# 添加标题和标签
plt.title("original and newdata after compensate(wz)")  # 设置整个图的标题
plt.xlabel("t(s)")  # x 轴标签
plt.ylabel("wz(deg/s)")  # y 轴标签
# 添加图例
plt.legend()
# 显示图形
plt.show()