from PIL import Image

# Paths
img1_path = "/home/vaibhav11/TE_plot/te_class_pie_grouped.png"
img2_path = "/home/vaibhav11/TE_plot/te_class_bar.png"
output_png = "/home/vaibhav11/TE_plot/combined2.png"
output_pdf = "/home/vaibhav11/TE_plot/combined2.pdf"

# Load images
img1 = Image.open(img1_path)
img2 = Image.open(img2_path)

# Optional: specify target height to resize both images proportionally
target_height = 600  # set this to desired height or None to keep original sizes

if target_height:
    def resize_proportional(img, target_h):
        w = int(img.width * target_h / img.height)
        return img.resize((w, target_h), resample=Image.LANCZOS)
    
    img1 = resize_proportional(img1, target_height)
    img2 = resize_proportional(img2, target_height)
else:
    # If no target height specified, resize both to the smaller height
    if img1.height != img2.height:
        new_height = min(img1.height, img2.height)
        img1 = img1.resize((int(img1.width * new_height / img1.height), new_height), resample=Image.LANCZOS)
        img2 = img2.resize((int(img2.width * new_height / img2.height), new_height), resample=Image.LANCZOS)

# Combine images side-by-side
combined_width = img1.width + img2.width
combined_image = Image.new("RGB", (combined_width, img1.height))
combined_image.paste(img1, (0, 0))
combined_image.paste(img2, (img1.width, 0))

# Save output
combined_image.save(output_png)
combined_image.save(output_pdf)

print(f"Combined image saved as:\n - {output_png}\n - {output_pdf}")

