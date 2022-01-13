from PIL import Image


def create_image(i, j):
  image = Image.new("RGB", (i, j), "white")
  return image


def get_pixel(image, i, j):
    # Inside image bounds?
    width, height = image.size
    if i > width or j > height:
      return None
    # Get Pixel
    pixel = image.getpixel((i, j))
    return pixel



#input:
nr_frames   = 10
im_w_cm     = 10.0
msk_open_cm = 0.1

#read:
im_ori      = Image.open('gif_project/hypnotic-shapes-moving-animated-gif-13.jpg')
im_w_nrpix  = im_ori.size[0] 
im_h_nrpix  = im_ori.size[1] 

#calc:
nrpix_w_percm   = im_w_nrpix/im_w_cm 
msk_open_nrpix  = int(msk_open_cm*nrpix_w_percm)
msk_clos_nrpix  = int((nr_frames - 1.0)*msk_open_nrpix)
msk_unit_nrpix  = int((nr_frames*msk_open_nrpix))
nr_msk_units    = int(im_w_nrpix/msk_unit_nrpix)


#resize picture to mek sure there is integer nr of msk units. that
#makes it all much better.
#make x-axis mask array: this is BEST AS WE CAN SHIFT IT FOR THE OTHE FRAMES!!!!


msk_1dw_arr = np.zeros(im_w_nrpix, dtype=np.int) 



print msk_open_nrpix, msk_clos_nrpix, nr_msk_units
#if we want, we can always change the initil res using the thumpnails routine...
print im_ori.format, im_ori.size, im_ori.mode


#input what frame it should be. start from 1.
frame_id = 1

#make a white copy of the image:
im_x    = create_image(im_w_nrpix, im_h_nrpix)
pixs_x  = im_x.load()

#pix i,j
#for j in range(0,im_h_nrpix)
        
    



im_ori.show()

#box = (100, 100, 400, 400)
#region = im.crop(box)

#print im.getpixel((10, 100))
#region.show()
























