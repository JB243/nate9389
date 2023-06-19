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

def counting_island(world: list)->int:
    class Node():
        def __init__(self, i, j): 
            self.i = i
            self.j = j
    
    class undirected_graph():     
        def __init__(self, V:list, E:list)->None:
            self.V = V[:]
            self.neighbor = {}
            for v in V:
                self.neighbor[v] = []
            for (v,w) in E:
                self.neighbor[v].append(w)

        def DFT_preorder(self)->int:
            count = 0
            if self.V:
                visited = {}
                for v in self.V:
                    visited[v]=False
                for v in self.V:
                    if not visited[v]: 
                        count += 1
                        self.__DFT__preorderHelp(visited, v)
            return count
                        
        def __DFT__preorderHelp(self, visited: list, v: int)->None:
            if not visited[v]:
                visited[v] = True
                for w in self.neighbor[v]:
                    self.__DFT__preorderHelp(visited, w)

    
    V = []
    E = []
    
    for i in range(len(world)):
        for j in range(len(world[0])):
            if world[i][j] == 1:
                V.append(Node(i, j))

    for v in V:
        for w in V:
            if w.i == v.i and w.j == v.j - 1:
                E.append((v,w))
            if w.i == v.i and w.j == v.j + 1:
                E.append((v,w))
            if w.i == v.i - 1 and w.j == v.j:
                E.append((v,w))
            if w.i == v.i + 1 and w.j == v.j:
                E.append((v,w))
    
    g = undirected_graph(V, E)
    
    return g.DFT_preorder()

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

def two_image_correlation_RG(img1_dir):

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

def two_image_correlation_RG_3D(my_dir):

    images = glob.glob(my_dir + "*.jpg")
    
    l_img1 = []
    l_img2 = []

    for k in range(len(images)):
        print(k)
        img1 = cv2.imread(images[k]) 
        img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
        for i in range(img1.shape[0]):
            for j in range(img1.shape[1]):
                if (img1[i:i+1, j:j+1, 0:1] + img1[i:i+1, j:j+1, 1:2] + img1[i:i+1, j:j+1, 2:3]) != 0:
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

def SSIM(x, y):
    # assumption : x and y are grayscale images with the same dimension

    import numpy as np
    
    def mean(img):
        return np.mean(img)
        
    def sigma(img):
        return np.std(img)
    
    def cov(img1, img2):
        img1_ = np.array(img1[:,:], dtype=np.float64)
        img2_ = np.array(img2[:,:], dtype=np.float64)
                        
        return np.mean(img1_ * img2_) - mean(img1) * mean(img2)
    
    K1 = 0.01
    K2 = 0.03
    L = 256 # when each pixel spans 0 to 255
   
    C1 = K1 * K1 * L * L
    C2 = K2 * K2 * L * L
    C3 = C2 / 2
        
    l = (2 * mean(x) * mean(y) + C1) / (mean(x)**2 + mean(y)**2 + C1)
    c = (2 * sigma(x) * sigma(y) + C2) / (sigma(x)**2 + sigma(y)**2 + C2)
    s = (cov(x, y) + C3) / (sigma(x) * sigma(y) + C3)
        
    return l * c * s

def mutual_information(img1, img2):
    import numpy as np
    import cv2
    import matplotlib.pyplot as plt
    
    # img1 and img2 are 3-channel color images
    
    a = img1[:,:,0:1].reshape(img1.shape[0], img1.shape[1])
    b = img2[:,:,0:1].reshape(img2.shape[0], img2.shape[1])
    
    hgram, x_edges, y_edges = np.histogram2d(
     a.ravel(),
     b.ravel(),
     bins=20
    )

    pxy = hgram / float(np.sum(hgram))
    px = np.sum(pxy, axis=1) # marginal for x over y
    py = np.sum(pxy, axis=0) # marginal for y over x
    px_py = px[:, None] * py[None, :] # Broadcast to multiply marginals

    nzs = pxy > 0 # Only non-zero pxy values contribute to the sum
    
    return np.sum(pxy[nzs] * np.log(pxy[nzs] / px_py[nzs]))

def spatial_featuremap(img, x, y, c):
    tsimg = np.zeros(img.shape[:2])    
    tsimg_row = y # np.array형 변수
    tsimg_col = x # np.array형 변수
    for rr, cc, t in zip(tsimg_row, tsimg_col, c):
        r, c = draw.circle(rr, cc, radius = 2.5) 
        tsimg[r,c]= t
    return tsimg

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

def binary_predict_by_1D_CNN (Xtr_rows, Ytr, Xte_rows, Yte):
    from tensorflow.keras.models import Sequential
    from tensorflow.keras.layers import Embedding, Dropout, Conv1D, GlobalMaxPooling1D, Dense
    from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint
    from tensorflow.keras.models import load_model
    
    vocab_size = 10000
    embedding_dim = 256 # 임베딩 벡터의 차원
    dropout_ratio = 0.3 # 드롭아웃 비율
    num_filters = 256 # 커널의 수
    kernel_size = 3 # 커널의 크기
    hidden_units = 128 # 뉴런의 수

    model = Sequential()
    model.add(Embedding(vocab_size, embedding_dim))
    model.add(Dropout(dropout_ratio))
    model.add(Conv1D(num_filters, kernel_size, padding='valid', activation='relu'))
    model.add(GlobalMaxPooling1D())
    model.add(Dense(hidden_units, activation='relu'))
    model.add(Dropout(dropout_ratio))
    model.add(Dense(1, activation='sigmoid'))
    
    model.summary()
    
    es = EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=3)
    mc = ModelCheckpoint('best_model.h5', monitor='val_acc', mode='max', verbose=1, save_best_only=True)

    model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['acc'])
    history = model.fit(Xtr_rows, Ytr, epochs=20, validation_data=(Xte_rows, Yte), callbacks=[es, mc])    

def binary_evaluate_by_1D_CNN(Xtr_rows, Ytr, Xte_rows, Yte):
    from tensorflow.keras.models import load_model
    binary_predict_by_1D_CNN(Xtr_rows, Ytr, Xte_rows, Yte)
    loaded_model = load_model('best_model.h5')
    print("\n 테스트 정확도: %.4f" % (loaded_model.evaluate(Xte_rows, Yte)[1]))

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

def GO(markers, max_terms = 6):
    # Reference
    ## https://github.com/mousepixels/sanbomics_scripts/blob/main/GO_in_python.ipynb
    ## https://snyk.io/advisor/python/goatools/functions/goatools.anno.genetogo_reader.Gene2GoReader
	
    # type(markers) = numpy.ndarray
    # dtype = object
    
    if len(markers) == 0:
        return
    
    from genes_ncbi_mus_musculus_proteincoding import GENEID2NT as GeneID2nt_mus
    from genes_ncbi_homo_sapiens_proteincoding import GENEID2NT as GeneID2nt_homo
    
    from goatools.base import download_go_basic_obo
    from goatools.base import download_ncbi_associations
    from goatools.obo_parser import GODag
    from goatools.anno.genetogo_reader import Gene2GoReader
    from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
    
    obo_fname = download_go_basic_obo()
    fin_gene2go = download_ncbi_associations()
    obodag = GODag("go-basic.obo")
    
    #run one time to initialize
    mapper = {}
    
    if markers[0] == markers[0].upper(): # Homo sapiens
        for key in GeneID2nt_homo:
            mapper[GeneID2nt_homo[key].Symbol] = GeneID2nt_homo[key].GeneID
    else:  # Mus musculus
        for key in GeneID2nt_mus:
            mapper[GeneID2nt_mus[key].Symbol] = GeneID2nt_mus[key].GeneID
        
    inv_map = {v: k for k, v in mapper.items()}
        
    if markers[0] == markers[0].upper():
        objanno = Gene2GoReader(fin_gene2go, taxids=[9606])
    else:
        objanno = Gene2GoReader(fin_gene2go, taxids=[10090])
    
    ns2assoc = objanno.get_ns2assc()

    goeaobj = ''
    if markers[0] == markers[0].upper():
        goeaobj = GOEnrichmentStudyNS(
            GeneID2nt_homo.keys(), # List of human protein-coding genes
            ns2assoc, # geneid/GO associations
            obodag, # Ontologies
            propagate_counts = False,
            alpha = 0.05, # default significance cut-off
            methods = ['fdr_bh']) # defult multipletest correction method
    else:
        goeaobj = GOEnrichmentStudyNS(
            GeneID2nt_mus.keys(), # List of mouse protein-coding genes
            ns2assoc, # geneid/GO associations
            obodag, # Ontologies
            propagate_counts = False,
            alpha = 0.05, # default significance cut-off
            methods = ['fdr_bh']) # defult multipletest correction method        
        
    GO_items = []

    temp = goeaobj.ns2objgoea['BP'].assoc
    for item in temp:
        GO_items += temp[item]
    
    temp = goeaobj.ns2objgoea['CC'].assoc
    for item in temp:
        GO_items += temp[item]
    
    temp = goeaobj.ns2objgoea['MF'].assoc
    for item in temp:
        GO_items += temp[item]
            
    def go_it(test_genes):
        print(f'input genes: {len(test_genes)}')
    
        mapped_genes = []
        for gene in test_genes:
            try:
                mapped_genes.append(mapper[gene])
            except:
                pass
        print(f'mapped genes: {len(mapped_genes)}')
        
        goea_results_all = goeaobj.run_study(mapped_genes)                
        goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]
                
        GO = pd.DataFrame(list(map(lambda x: [x.GO, x.goterm.name, x.goterm.namespace, x.p_uncorrected, x.p_fdr_bh,\
                       x.ratio_in_study[0], x.ratio_in_study[1], GO_items.count(x.GO), list(map(lambda y: inv_map[y], x.study_items)),\
                       ], goea_results_sig)), columns = ['GO', 'term', 'class', 'p', 'p_corr', 'n_genes',\
                                                        'n_study', 'n_go', 'study_genes'])

        GO = GO[GO.n_genes > 1]        
        return GO
        
    df = go_it(markers)    
    df['per'] = df.n_genes/df.n_go
    
    if df.shape[0] > 0:
        df = df.sort_values(by=['p_corr'], axis=0, ascending=True)
        
        if df.shape[0] > max_terms:
            df = df[0:max_terms]

        import matplotlib.pyplot as plt
        import matplotlib as mpl
        from matplotlib import cm
        import seaborn as sns
        import textwrap
    
        fig, ax = plt.subplots(figsize = (0.5, 2.75))
        cmap = mpl.cm.bwr_r
        norm = mpl.colors.Normalize(vmin = df.p_corr.min(), vmax = df.p_corr.max())
        mapper = cm.ScalarMappable(norm = norm, cmap = cm.bwr_r)
        cbl = mpl.colorbar.ColorbarBase(ax, cmap = cmap, norm = norm, orientation = 'vertical')
    
        plt.figure(figsize = (2,4))
        ax = sns.barplot(data = df, x = 'per', y = 'term', palette = mapper.to_rgba(df.p_corr.values))
        ax.set_yticklabels([textwrap.fill(e, 22) for e in df['term']])
        plt.show()

def GO_(gene_list):
    # reference : https://gseapy.readthedocs.io/en/latest/gseapy_example.html
    
    import gseapy
    from gseapy import barplot, dotplot
    
    if gene_list[0] == gene_list[0].upper():
        enr = gseapy.enrichr(gene_list=gene_list, # or "./tests/data/gene_list.txt",
                     gene_sets=['GO_Biological_Process_2018','GO_Cellular_Component_2018', 'GO_Molecular_Function_2018'],
                     organism='human', # don't forget to set organism to the one you desired! e.g. Yeast
                     outdir=None, # don't write to disk
                    )
    else:
        enr = gseapy.enrichr(gene_list=gene_list, # or "./tests/data/gene_list.txt",
             gene_sets=['GO_Biological_Process_2018','GO_Cellular_Component_2018', 'GO_Molecular_Function_2018'],
             organism='mouse', # don't forget to set organism to the one you desired! e.g. Yeast
             outdir=None, # don't write to disk
            )

    try:
        ax = dotplot(enr.results,
                      column="Adjusted P-value",
                      x='Gene_set', # set x axis, so you could do a multi-sample/library comparsion
                      size=10,
                      top_term=5,
                      figsize=(3,5),
                      title = "GO plot",
                      xticklabels_rot=45, # rotate xtick labels
                      show_ring=True, # set to False to revmove outer ring
                      marker='o',
                     )
    except:
        print("ValueError: Warning: No enrich terms when cutoff = 0.05")

def ensembl_to_gene(ensembl_list):
    ar = []

    human = pd.read_csv("https://blog.kakaocdn.net/dn/29YTj/btrS5iG9QOH/Di6RQKxHOPDii7EjkdHN30/human_genes_36601.tsv?attach=1&knm=tfile.tsv", sep = '\t', header = None)
    mouse = pd.read_csv("https://blog.kakaocdn.net/dn/wkjwJ/btrS1QSgrpD/VS8ELANCQyeZAA3vL8JQP0/mouse_genes_32285.tsv?attach=1&knm=tfile.tsv", sep = '\t', header = None)

    for i in range(len(ensembl_list)):
        if 'ENSG' in ensembl_list[i]:  # human gene
            index = human[0].tolist().index(ensembl_list[i])
            ar.append(human[1][index])
        elif 'ENSMUSG' in ensembl_list[i]:  # mouse gene
            index = mouse[0].tolist().index(ensembl_list[i])
            ar.append(mouse[1][index])

    return(ar)

def gene_to_ensembl(gene_list):
    ar = []

    human = pd.read_csv("https://blog.kakaocdn.net/dn/29YTj/btrS5iG9QOH/Di6RQKxHOPDii7EjkdHN30/human_genes_36601.tsv?attach=1&knm=tfile.tsv", sep = '\t', header = None)
    mouse = pd.read_csv("https://blog.kakaocdn.net/dn/wkjwJ/btrS1QSgrpD/VS8ELANCQyeZAA3vL8JQP0/mouse_genes_32285.tsv?attach=1&knm=tfile.tsv", sep = '\t', header = None)

    for i in range(len(gene_list)):
        if gene_list[i] == gene_list[i].upper():  # human gene
            index = human[1].tolist().index(gene_list[i])
            ar.append(human[0][index])
        else:  # mouse gene
            index = mouse[1].tolist().index(gene_list[i])
            ar.append(mouse[0][index])

    return(ar)

def find_idx(my_list, e):
    for i in range(len(my_list)):
        if my_list[i] == e:
            return i
    return -1

def human_to_mouse(human_gene:list):
    hom = pd.read_csv("https://blog.kakaocdn.net/dn/cVeqsA/btrS1JMnxyX/HtVhPmqtxdgt7LQlGkeql0/HOM_MouseHumanSequence.csv?attach=1&knm=tfile.csv")

    mouse_gene = []

    for i in range(len(human_gene)):
        index = find_idx(hom['Symbol'], human_gene[i])
        DB_Class_Key = hom['DB Class Key'][index]
        hom_ = hom[hom['DB Class Key'] == DB_Class_Key]
        hom__ = hom_[hom_['Common Organism Name'] == 'mouse, laboratory']
        mouse_gene.append(np.array(hom__['Symbol'])[0])

    return mouse_gene

def mouse_to_human(mouse_gene:list):
    hom <- read.csv("https://blog.kakaocdn.net/dn/cVeqsA/btrS1JMnxyX/HtVhPmqtxdgt7LQlGkeql0/HOM_MouseHumanSequence.csv?attach=1&knm=tfile.csv")

    human_gene = []

    for i in range(len(mouse_gene)):
        index = find_idx(hom['Symbol'], mouse_gene[i])
        DB_Class_Key = hom['DB Class Key'][index]
        hom_ = hom[hom['DB Class Key'] == DB_Class_Key]
        hom__ = hom_[hom_['Common Organism Name'] == 'human']
        mouse_gene.append(np.array(hom__['Symbol'])[0])

    return human_gene
