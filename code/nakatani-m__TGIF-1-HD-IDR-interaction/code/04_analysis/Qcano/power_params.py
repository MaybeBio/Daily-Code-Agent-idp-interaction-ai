import numpy as np
from scipy.ndimage import uniform_filter

def rw_params(path,power=1):
    def smooth(values,coords,dims):
        #配列を画像へ変換
        arr1,arr_original = np.zeros(dims),np.zeros(dims)
        for c, v in zip(coords[:-1], values[:-1]):
            arr_original[c[0]-1,c[1]-1, c[2]-1] = v
            arr1[c[0]-1,c[1]-1,c[2]-1] = 1
        
        #スムージング処理
        arr1_smooth = uniform_filter(arr1, size=3, mode='constant', cval=0)
        arr_smooth = uniform_filter(arr_original, size=3, mode='constant', cval=0)
        with np.errstate(divide='ignore', invalid='ignore'):
            arr_smooth = (arr_smooth*27 - arr_original)/(27*arr1_smooth - 1)
        
        #画像を配列へ変換
        tmp = []
        for c in coords[:-1]:
            tmp.append(arr_smooth[c[0]-1,c[1]-1,c[2]-1])
        tmp.append(values[-1])
        return np.array(tmp)
    
    #ファイルの読み込み
    with open(path+"params_next.dat") as f:
        lines = f.readlines()
    dim = int(lines[1])
    dims = []
    tmp_n = 2
    for _ in range(dim):
        size = int(lines[tmp_n])
        dims.append(size)
        tmp_n += size + 1
    values,coords,tmp = [],[],[]
    for line in lines:
        txt = line.split()
        if len(txt) >= dim + 1:
            coords.append([int(txt[i]) for i in range(dim)])
            values.append(float(txt[dim]))
        elif txt[0]!="END":
            tmp.append(line)
    values = np.array(values)

    #データの比較
    values_smooth = smooth(values,coords,dims)
    values_result = np.minimum(values_smooth*100,values**power)
    
    #データの書き出し
    with open("params_next_power.dat", "w") as fo:
        fo.writelines(tmp)
        for c, v in zip(coords[:-1], values_result[:-1]):
            fo.write(f"{c[0]} {c[1]} {c[2]} {v}\n")
        fo.write("{} {} {} {}\n".format(coords[-1][0],coords[-1][1], coords[-1][2], values[-1]))
        fo.write("END")

    #修正データの保存
    with open("param_changes.log", "w") as fo:
        for c, v, v_r in zip(coords[:-1], values[:-1],values_result[:-1]):
            if v != v_r:
                fo.write(f"{c[0]} {c[1]} {c[2]} {v}\n")   

def main():
    path="./"
    rw_params(path)

if __name__=="__main__":
    main()
