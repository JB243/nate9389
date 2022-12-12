from keras.layers import DepthwiseConv2D
from keras.layers import Input
from keras.models import Model
from get_model import *

def squared_error(y_true, y_pred):
    return K.sum(K.square(y_pred - y_true))

def define_denoising_model(height, width, channel=1):
    num_u = [32, 32, 32, 32, 32]#original [128, 128, 128, 128, 128]
    num_d = [32, 32, 32, 32, 32] #[128, 128, 128, 128, 128]
    kernel_u = [3, 3, 3, 3, 3]
    kernel_d = [3, 3, 3, 3, 3]
    num_s = [4, 4, 4, 4, 4]
    kernel_s = [1, 1, 1, 1, 1]
    lr = 0.01
    inter = "bilinear"

    model = define_model(num_u, num_d, kernel_u, kernel_d, num_s, kernel_s, height, width, inter, lr,
                         input_channel=channel)
    model.compile(loss=squared_error, optimizer=Adam(lr=lr))

    return model


def denoising(image, init_image=None, num_iter = 500, verbose=0):
    height, width = image.shape[:2]
    if init_image.all()==None:
        channels = 1
        model = define_denoising_model(height, width, channel=channels)
        input_noise = np.random.uniform(0, 0.1, (1, height, width, 1))
    else:
        channels = init_image.shape[2]
        model = define_denoising_model(height, width, channel=channels)
        input_noise = init_image[None, :, :, :]
        
    if verbose:
        plt.ion()
    for i in range(num_iter):
        model.train_on_batch(input_noise + np.random.normal(0, 1/30.0, (height, width, channels)), image[None, :, :, :])
        if verbose:
            if i % 100==0:
                outimg = model.predict(input_noise)[0]
                plt.imshow(outimg[:,:,0])
                plt.title(("Iteration:"+str(i+1)))
                plt.show()
    return model.predict(input_noise)[0]
    
    #return np.clip(model.predict(input_noise)[0], 0, 255).astype(np.uint8)