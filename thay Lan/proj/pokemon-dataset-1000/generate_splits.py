import os
import shutil
import random

# Paths
base_path = "."
dataset_path = os.path.join(base_path, "dataset")
train_path = os.path.join(base_path, "train")
val_path = os.path.join(base_path, "val")
test_path = os.path.join(base_path, "test")

# Create directories for train, val, and test splits
os.makedirs(train_path, exist_ok=True)
os.makedirs(val_path, exist_ok=True)
os.makedirs(test_path, exist_ok=True)

# Train-Val-Test split ratios
train_ratio = 0.8
val_ratio = 0.1
test_ratio = 0.1

# Iterate through class folders
for class_name in os.listdir(dataset_path):
    class_folder = os.path.join(dataset_path, class_name)
    if os.path.isdir(class_folder):  # Check if it's a directory
        # Get all image files in the class folder
        images = [os.path.join(class_folder, img) for img in os.listdir(class_folder) if img.endswith('.png')]

        # Shuffle images for randomness
        random.shuffle(images)

        # Calculate split sizes
        train_size = int(len(images) * train_ratio)
        val_size = int(len(images) * val_ratio)

        # Split images
        train_images = images[:train_size]
        val_images = images[train_size:train_size + val_size]
        test_images = images[train_size + val_size:]

        # Create class subfolders in train, val, and test directories
        os.makedirs(os.path.join(train_path, class_name), exist_ok=True)
        os.makedirs(os.path.join(val_path, class_name), exist_ok=True)
        os.makedirs(os.path.join(test_path, class_name), exist_ok=True)

        # Move images to respective directories
        for img in train_images:
            shutil.copy(img, os.path.join(train_path, class_name))
        for img in val_images:
            shutil.copy(img, os.path.join(val_path, class_name))
        for img in test_images:
            shutil.copy(img, os.path.join(test_path, class_name))

        print(f"Processed class: {class_name}")

print("Dataset split into train, val, and test folders.")
