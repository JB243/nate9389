def SUM(a, b):
    return a+b

def is_contain(given_str, partial_str):    
    s = ''
    for char in given_str:
        s = s + char
        if s.endswith(partial_str):
            return True
    
    return False

def str_omit_forward(_str, n):
    s = ''
    count = 0
    for char in _str:
        count = count + 1
        if count > n:
            s = s + char
    return s

def str_omit_backward(_str, n):
    s = ''
    count = 0
    for char in _str:
        count = count + 1
        if count <= len(_str) - n:
            s = s + char
    return s

def is_element_in_list(e, l):
    for i in range(len(l)):
        if e == l[i]:
            return True
    return False

def is_element_in_list(e, l):
    for i in range(len(l)):
        if e == l[i]:
            return True
    return False

def is_element_in_list(e, l):
    for i in range(len(l)):
        if e == l[i]:
            return True
    return False

def read_txt(txt_dir):

    import csv
    import sys
    
    l = []
    f = open(txt_dir, 'r', encoding = 'utf-8')
 
    for line in csv.reader(f):
        ### if you want elements of list ###
        l.append(line) 
		
        ### if you want elements of str ###
        # l.append(''.join(line)) 
    
    f.close()
    return l

def two_image_correlation_RG(img1_dir):

    # img1 and img2 should be same in size
    img1 = cv2.imread(img1_dir) # RGB image
    img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)

    l_img1 = []
    l_img2 = []

    for i in range(img1.shape[0]):
        for j in range(img1.shape[1]):
            l_img1.append((
                img1[i:i+1, j:j+1, 0:1] 
            ).item())

            l_img2.append((
                img1[i:i+1, j:j+1, 1:2] 
            ).item())

    print("brightness of img1")
    print(np.mean(l_img1))
    print("brightness of img2")
    print(np.mean(l_img2))

    print("img1-img2 correlation")
    print(scipy.stats.pearsonr(l_img1, l_img2) )

    plt.scatter(l_img1, l_img2, alpha = 0.01)
    plt.xlabel('brightness of img1')
    plt.ylabel('brightness of img2')

def two_image_correlation(img1_dir, img2_dir):

    import cv2
    import matplotlib.pyplot as plt
    from skimage.util import view_as_windows
    import numpy as np
    import scipy.stats
    
    # img1 and img2 should be same in size
    img1 = cv2.imread(img1_dir) # RGB image
    img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
    img2 = cv2.imread(img2_dir) # RGB image
    img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)
		
    l_img1 = []
    l_img2 = []
		
    for i in range(img1.shape[0]):
        for j in range(img1.shape[1]):
            l_img1.append((
                0.2989*img1[i:i+1, j:j+1, 0:1] + 
                0.5870*img1[i:i+1, j:j+1, 1:2] + 
                0.1140*img1[i:i+1, j:j+1, 2:3]
            ).item())
		        
            l_img2.append((
                0.2989*img2[i:i+1, j:j+1, 0:1] + 
                0.5870*img2[i:i+1, j:j+1, 1:2] + 
                0.1140*img2[i:i+1, j:j+1, 2:3]
            ).item())
		
    print("brightness of img1")
    print(np.mean(l_img1))
    print("brightness of img2")
    print(np.mean(l_img2))
		
    print("img1-img2 correlation")
    print(scipy.stats.pearsonr(l_img1, l_img2) )
		
    plt.scatter(l_img1, l_img2, alpha = 0.05)
    plt.xlabel('brightness of img1')
    plt.ylabel('brightness of img2')
    # sns.regplot(l_b_CD31, l_b_Lipo, alpha = 0.05)
    
def RGBtoGray(img):
    import numpy as np

    out = np.empty((img.shape[0], img.shape[1]))

    for i in range(img.shape[0]):
        for j in range(img.shape[1]):
            out[i:i+1, j:j+1] = ( 0.2989*img[i:i+1, j:j+1, 0:1] + 
                                  0.5870*img[i:i+1, j:j+1, 1:2] + 
                                  0.1140*img[i:i+1, j:j+1, 2:3]
                                ).item()
      
    return out
    
def recolor(img, pre, post): 
    for i in range(img.shape[0]):
        for j in range(img.shape[1]):
            if img[i:i+1, j:j+1] == pre:
                img[i:i+1, j:j+1] = post
    return(img)

def KNN_predict(Xtr_rows, Ytr, Xte_rows, dist_metric='L2'):
    import numpy as np
    
    class NearestNeighbor(object):
        def __init__(self):
            pass

        def train(self, X, Y):
            self.Xtr = X
            self.Ytr = Y

        def predict(self, X, dist_metric=dist_metric):
            num_test = X.shape[0]
            Ypred = np.zeros(num_test, dtype = self.Ytr.dtype)
        
            for i in range(num_test):
                if dist_metric=='L1': ## L1 distance
                    distances = np.sum(np.abs(self.Xtr - X[i, :]), axis = 1) 
                elif dist_metric=='L2': ## L2 distance
                    distances = np.sum(np.square(self.Xtr - X[i, :]), axis = 1) 
                elif dist_metric=='dot': ## dot distance
                    distances = np.dot(self.Xtr, X[i, :])
            
                min_index = np.argmin(distances)
                Ypred[i] = self.Ytr[min_index]
            
                if i%100 == 0:
                    print(i)
            return Ypred
        
    nn = NearestNeighbor()
    nn.train(Xtr_rows, Ytr)
    
    Yte_predict = nn.predict(Xte_rows, dist_metric)
    return Yte_predict


def KNN_evaluate(Xtr_rows, Ytr, Xte_rows, Yte, dist_metric='L2'):
    Yte_predict = KNN_predict(Xtr_rows, Ytr, Xte_rows, dist_metric)
    print ('accuracy: %f' + str(np.mean(Yte_predict == Yte)) )
    # L1 : accuracy = 47.29%
    # L2 : accuracy = 27.24%
    # dot : accuracy = 9.9%

def features_512D_from_image_by_CNN (img):
    # Image Patch
    img_re = cv2.resize(img, dsize = (32, 32))
    tspatches = []
    tspatches.append(img_re)
    tspatches.append(img_re) # intentional duplicate
    tspatches = np.asarray(tspatches)    
    
    # Pretrained model
    pretrained_model = vgg16.VGG16(weights='imagenet', include_top = False, pooling='avg', input_shape = (32,32,3))
    X_in = tspatches.copy()
    X_in = vgg16.preprocess_input(X_in)
    pretrained_model.trainable = False
    pretrained_model.summary()
    
    ts_features = pretrained_model.predict(X_in)
    
    return ts_features[0]
