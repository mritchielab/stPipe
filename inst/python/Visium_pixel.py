import cv2
import numpy as np
import pandas as pd
import argparse
import os

# setup input and resize
def read_img(img_path, resized_shape=(10000, 10000)):
    # read image
    org_img = cv2.imread(img_path)
    if org_img is None:
        raise FileNotFoundError(f"Image file not found or cannot be read: {img_path}")
    
    # resize but keep the resize scale for later use
    resize_scale = (resized_shape[1] / org_img.shape[1], resized_shape[0] / org_img.shape[0])
    img = cv2.resize(org_img, (resized_shape[0], resized_shape[1]))
    # Convert to grayscale (fiducials are often distinct in grayscale)
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    # Apply Gaussian blur to reduce noise
    blurred = cv2.GaussianBlur(gray, (9, 9), 2)
    # Threshold the image to get a binary image with gaussian adaptive thresholding
    thres_param = round(resized_shape[0] / 400)
    if thres_param % 2 == 0:
        thres_param += 1
    thresholded = cv2.adaptiveThreshold(blurred, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY, thres_param, 1.5)
    return thresholded, resize_scale, org_img

def process_query(img_fn, ref_p, ref_spots, fiducials, plot, output_mapping_image):
    query_img, resize_scale, org_img = read_img(img_fn)
    query = np.float32([find_fiducial(query_img, fid) for fid in fiducials])
    matrix = cv2.getPerspectiveTransform(ref_p, query)
    transformed_points = cv2.perspectiveTransform(ref_spots, matrix).reshape(-1, 2)
    transformed_points = transformed_points / resize_scale
    transformed_points = np.round(transformed_points).astype(int)
    
    # Draw transformed points on the original image
    for point in transformed_points:
        cv2.circle(org_img, tuple(point), org_img.shape[0] // 200, (0, 255, 0), -1)

    cv2.imwrite(output_mapping_image, org_img)
    
    return transformed_points, resize_scale

def find_fiducial(image, fiducial):
    result = cv2.matchTemplate(image, fiducial, cv2.TM_CCOEFF_NORMED)
    _, _, _, max_loc = cv2.minMaxLoc(result)
    return max_loc

def resize_coord(pos, org_w, org_h, org_scale, res_w, res_h):
    x, y = pos
    x_new = int(x * res_w / org_w * org_scale)
    y_new = int(y * res_h / org_h * org_scale)
    return x_new, y_new

def main(test_img_fn, output_dir):
    package_dir = os.getenv("R_PACKAGE_DIR", default=".")
    ref_img_fn = os.path.join(package_dir, "extdata/tissue_hires_image.png")
    ref_pos_fn = os.path.join(package_dir, "extdata/ref_pos.csv")
    ref_scale = 0.06300403

    # load reference
    ref_img, _, _ = read_img(ref_img_fn)

    # pre-define reference area
    tmpscale = ref_img.shape[0] / 2000
    y_range1 = slice(int(80 * tmpscale), int(260 * tmpscale))
    y_range2 = slice(int(1650 * tmpscale), int(1850 * tmpscale))
    x_range1 = slice(int(104 * tmpscale), int(300 * tmpscale))
    x_range2 = slice(int(1700 * tmpscale), int(1900 * tmpscale))

    fiducials = [
        ref_img[y_range1, x_range1].copy(),
        ref_img[y_range2, x_range1].copy(),
        ref_img[y_range1, x_range2].copy(),
        ref_img[y_range2, x_range2].copy()
    ]

    ref_p = np.float32([
        [y_range1.start, x_range1.start],
        [y_range2.start, x_range1.start],
        [y_range1.start, x_range2.start],
        [y_range2.start, x_range2.start]
    ])

    org_ref_img = cv2.imread(ref_img_fn)
    original_width, original_height = org_ref_img.shape[1], org_ref_img.shape[0]
    resized_width, resized_height = ref_img.shape[1], ref_img.shape[0]

    # read the reference spot coordinates
    ref_spots = pd.read_csv(ref_pos_fn).iloc[:, -2:].values.astype(np.float32)
    ref_spots = np.array([resize_coord(x, original_width, original_height, ref_scale, resized_width, resized_height)
                          for x in ref_spots], dtype='float32')
    ref_spots = np.array([ref_spots], dtype='float32')

    # Prepare output directory
    os.makedirs(output_dir, exist_ok=True)

    # Set the path for the mapping image
    output_mapping_image_path = os.path.join(output_dir, "mapping_output.jpg")
    
    # Process the query and generate results
    out_p_test, resize_scale = process_query(test_img_fn, ref_p, ref_spots, fiducials, plot=True, output_mapping_image=output_mapping_image_path)

    # Save result as CSV file
    output_csv_path = os.path.join(output_dir, "mapped_pixel.csv")
    pd.DataFrame(out_p_test, columns=['y', 'x']).to_csv(output_csv_path, index=False)
    
    print(f"Computed pixel results saved to {output_csv_path}")
    print(f"Mapping image saved to {output_mapping_image_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Image Processing Script for 10X Visium")
    parser.add_argument("test_img_fn", type=str, help="Path to the tiff image file")
    parser.add_argument("output_dir", type=str, help="Path to the directory to save all output results (pixel in CSV file format and mapping image in png file format)")
    args = parser.parse_args()
    main(args.test_img_fn, args.output_dir)
