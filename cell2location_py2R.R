library(Seurat)
library(dplyr)
library(stringr) # Used for easier string manipulation

# ================= 1. Read and fix column names =================
# Read CSV (assuming you have already read it in)
c2l_data <- read.csv("D:/Data_us/spatial_dataus_out/cell2location_results_for_R.csv", row.names = 1)

# Backup the original data just in case
c2l_clean <- c2l_data

# --- Core: Logic to fix column names ---
# R's logic is usually: "CD8+ Trm" -> "CD8..Trm" (two dots)
#                "Treg Th1" -> "Treg.Th1" (one dot)

original_colnames <- colnames(c2l_clean)
cols_to_remove <- c("NA.", "NA", "Unknown", "unassigned") # You can put everything you don't want here
c2l_clean <- c2l_clean[, !colnames(c2l_clean) %in% cols_to_remove]
# Step 1: Replace two consecutive dots ".." back to "+ " (plus and space)
# Note: Need to escape the dot
new_colnames <- gsub("\\.\\.", "+ ", original_colnames)

# Step 2: Replace remaining single dots "." back to " " (space)
new_colnames <- gsub("\\.", " ", new_colnames)

# Step 3: Apply new column names
colnames(c2l_clean) <- new_colnames

# Check the correction effect
print(">>> Column name repair example:")
print(colnames(c2l_clean)) # Check if the first few have changed back

# Continue modifying
rename_map <- c(
  "CD4+ CTL+ GNLY+ "         = "CD4+ CTL (GNLY+)",
  "CD4+ T_ISG+ STAT1+ "      = "CD4+ T_ISG (STAT1+)",
  "CD4+ Tem+ Cytotoxic "     = "CD4+ Tem (Cytotoxic)",
  "CD4+ Tfh+ CXCL13+ "       = "CD4+ Tfh (CXCL13+)",
  "CD8+ Pre Trm"             = "CD8+ Pre-Trm",             # Fixed hyphen
  "CD8+ Tcm+ Stem like "     = "CD8+ Tcm (Stem-like)",
  "CD8+ Tem+ GZMB+ "         = "CD8+ Tem (GZMB+)",
  "CD8+ Tex+ CXCL13+ "       = "CD8+ Tex (CXCL13+)",
  "CD8+ Tex+ GZMK+ High TOX " = "CD8+ Tex (GZMK+ High TOX)",
  "CD8+ Trm+ Resting "       = "CD8+ Trm (Resting)",
  "Epithelial Keratinocytes" = "Epithelial/Keratinocytes", # Fixed slash
  "Gd T"                     = "Gd-T",                     # Fixed hyphen
  "Melanoma+ HLA DRA+ "      = "Melanoma (HLA-DRA+)",
  "Naive Stem like Program"  = "Naive/Stem-like Program"   # Fixed slash and hyphen
)

# 2. Execute replacement
# Check column names of c2l_clean, if they appear in the dictionary, replace them with new names
current_names <- colnames(c2l_clean)

for (bad_name in names(rename_map)) {
  # If this garbled name exists in the data
  if (bad_name %in% current_names) {
    # Find its position and rename
    colnames(c2l_clean)[colnames(c2l_clean) == bad_name] <- rename_map[[bad_name]]
  }
}

# 3. Check the first few column names again to confirm columns that do not need fixing (like B Cells) remain unchanged, and those needing fixes are reverted to original
print(">>> Column name repair complete, check as follows:")
print(head(colnames(c2l_clean), 20)) # See if there are parentheses now

rm(t_cells_clean)

# ================= 2. Split data and clean Barcode =================

# Split the large table into a list by the sample column
df_list <- split(c2l_clean, c2l_clean$sample)

# Define the mapping relationship between sample names and Seurat object names
# Key is the sample name in CSV, Value is the object name in your R environment
obj_mapping <- list(
  "spatial_T13" = "test_data_T13",
  "spatial_T5"  = "test_data_T5",
  "spatial_T14" = "test_data_T14"
)

# ================= 3. Loop processing: Clean Barcode and inject into Seurat =================

# This step will automatically update test_data_T13, test_data_T5, test_data_T14 in your environment
# We don't need to write them manually one by one

for (sample_id in names(df_list)) {
  
  # 1. Get the corresponding Seurat object name
  seurat_obj_name <- obj_mapping[[sample_id]]
  
  if (exists(seurat_obj_name)) {
    message(paste(">>> Processing sample:", sample_id, "-> Corresponding object:", seurat_obj_name))
    
    # Get the C2L data for this sample
    current_df <- df_list[[sample_id]]
    
    # 2. Clean Barcode (Row names)
    # Logic: Remove "-sample_id" at the end
    # E.g., "AAAC...-1-spatial_T13" -> "AAAC...-1"
    # Use regex: replace the last hyphen and everything after it
    # Note: If your original barcode doesn't have -1 itself, this might need adjustment.
    # Usually Seurat's Barcode format is Sequence-1.
    
    # Here we precisely remove the "-sample_id" segment
    suffix_pattern <- paste0("-", sample_id, "$") 
    clean_barcodes <- gsub(suffix_pattern, "", rownames(current_df))
    rownames(current_df) <- clean_barcodes
    
    # 3. Prepare expression matrix (Features x Cells)
    # Remove metadata columns (x, y, sample), keep only cell abundance
    meta_cols <- c("x", "y", "sample")
    # Find columns that are not metadata
    abundance_data <- current_df[, !colnames(current_df) %in% meta_cols]
    
    # Seurat requires the matrix to be: rows=genes (cell types), columns=cells (Barcodes)
    # So transposition t() is needed
    c2l_matrix <- t(as.matrix(abundance_data))
    
    # 4. Get target Seurat object
    target_obj <- get(seurat_obj_name)
    
    # 5. Safety check: take intersection to prevent errors
    # Only cells present in both Seurat and C2L results will be added
    common_cells <- intersect(colnames(target_obj), colnames(c2l_matrix))
    
    if (length(common_cells) > 0) {
      # Subset matrix
      c2l_matrix_subset <- c2l_matrix[, common_cells, drop = FALSE]
      
      # 6. Create new Assay and add
      # We name this Assay "cell2location"
      c2l_assay <- CreateAssayObject(data = c2l_matrix_subset)
      target_obj[["cell2location"]] <- c2l_assay
      
      # (Optional) Set as default Assay for immediate plotting, comment out if not wanted
      # DefaultAssay(target_obj) <- "cell2location"
      
      # 7. Write the updated object back to the global environment
      assign(seurat_obj_name, target_obj)
      
      message(paste("    ✅ Successfully added annotation info of", nrow(c2l_matrix_subset), "cell types to", length(common_cells), "Spots."))
      
    } else {
      warning(paste("    ❌ Warning: No matching Barcode found in", seurat_obj_name, "! Please check Barcode format."))
      # Print the barcodes on both sides to facilitate debugging
      print("    C2L Barcode example:")
      print(head(colnames(c2l_matrix)))
      print("    Seurat Barcode example:")
      print(head(colnames(target_obj)))
    }
    
  } else {
    warning(paste("    ⚠️ Seurat object not found in environment:", seurat_obj_name))
  }
}

# ================= 4. Verify results =================

message("\n>>> Verifying if test_data_T13 contains cell2location Assay:")
if(exists("test_data_T13")) {
  print(Assays(test_data_T13))
  # Try printing the first few rows of data
  # print(GetAssayData(test_data_T13, assay = "cell2location")[1:3, 1:3])
}


# Plotting
library(Seurat)
library(ggplot2)
library(viridis)
library(patchwork)

# 1. Prepare sample list
# Ensure these objects are already in your environment
obj_list <- list(
  "T13" = test_data_T13,
  "T5"  = test_data_T5,
  "T14" = test_data_T14
)

# 2. Define the cell types you want to plot (must exactly match the repaired names)
target_cell_1 <- "CD8+ Trm"
target_cell_2 <- "Treg Th1 Program"

# ================= Plotting loop =================
for (name in names(obj_list)) {
  obj <- obj_list[[name]]
  
  # Switch default Assay to cell2location
  DefaultAssay(obj) <- "cell2location"
  
  message(paste(">>> Plotting sample:", name))
  
  # Plot CD8+ Trm
  p1 <- SpatialFeaturePlot(obj, features = target_cell_1, 
                           pt.size.factor = 1.6, # If dots are too small, increase this (e.g., 2.5 or 3.0)
                           alpha = c(0.1, 1),    # Transparency range
                           stroke = 0) +         # Remove point borders
    scale_fill_viridis_c(option = "magma") +     # Use magma colors (black-red-bright yellow)
    ggtitle(paste(name, target_cell_1))
  
  # Plot Treg
  p2 <- SpatialFeaturePlot(obj, features = target_cell_2, 
                           pt.size.factor = 1.6, 
                           alpha = c(0.1, 1), 
                           stroke = 0) +
    scale_fill_viridis_c(option = "plasma") +
    ggtitle(paste(name, target_cell_2))
  
  # Display combined plots
  print(p1 | p2)
}
# Change colors
library(scales) # Need to load this package to handle color gradients

# ================= Plotting loop (High contrast version) =================
for (name in names(obj_list)) {
  obj <- obj_list[[name]]
  DefaultAssay(obj) <- "cell2location"
  
  message(paste(">>> Plotting sample:", name))
  
  # --- 1. Plot CD8+ Trm (Using bright yellow-orange-red gradient) ---
  # Yellow is most conspicuous on a purple background
  p1 <- SpatialFeaturePlot(obj, features = target_cell_1, 
                           pt.size.factor = 2.8, # Increase point size (1.6 might have been too small before)
                           alpha = c(0, 0.9),    # [Key] Low expression set to 0 (fully transparent), high expression 0.9
                           stroke = 0) +
    # Use custom gradient: transparent/light grey -> bright yellow -> bright red
    scale_fill_gradientn(colours = c(alpha("grey90", 0.1), "yellow", "red"),
                         values = c(0, 0.4, 1)) + 
    ggtitle(paste(name, "CD8+ Trm")) +
    theme(legend.position = "right", 
          legend.text = element_text(size=10),
          plot.title = element_text(size=14, face="bold"))
  
  # --- 2. Plot Treg (Using bright cyan-blue gradient) ---
  # Cyan green forms a complementary color with purplish red background, very clear
  p2 <- SpatialFeaturePlot(obj, features = target_cell_2, 
                           pt.size.factor = 2.8, 
                           alpha = c(0, 0.9), 
                           stroke = 0) +
    # Use custom gradient: transparent/light grey -> bright cyan -> dark blue
    scale_fill_gradientn(colours = c(alpha("grey90", 0.1), "cyan", "blue"),
                         values = c(0, 0.4, 1)) +
    ggtitle(paste(name, "Th1 Program Treg ")) +
    theme(legend.position = "right",
          legend.text = element_text(size=10),
          plot.title = element_text(size=14, face="bold"))
  
  # Display combined plots
  print(p1 | p2)
}
# Advanced color scheme
library(scales)

# ================= 💎 Option A: Premium Gold & Glacier Blue 💎 =================

# 1. Define Gold series for Trm (Golden/Amber)
# Also from transparent -> bright gold -> deep orange, retain brightness but remove harsh pure yellow
col_trm_low  <- alpha("#FFF7BC", 0.1) # Extremely light beige (base)
col_trm_mid  <- "#FFC125"             # Royal Gold (Goldenrod1) - main highlight color
col_trm_high <- "#FF4500"             # Orange red (OrangeRed1) - finish with this to add depth

# 2. Define Glacier series for Treg (Glacial/Teal)
# From transparent -> Tiffany blue -> Sapphire blue
col_treg_low  <- alpha("#E0F3DB", 0.1) # Extremely light mint (base)
col_treg_mid  <- "#00E5EE"             # Turquoise (Turquoise2) - highly translucent
col_treg_high <- "#00688B"             # Deep peacock blue - adds sense of depth

# ================= Plotting loop =================
for (name in names(obj_list)) {
  obj <- obj_list[[name]]
  DefaultAssay(obj) <- "cell2location"
  
  message(paste(">>> Plotting sample (Premium texture version):", name))
  
  # --- 1. Plot CD8+ Trm (Liquid gold texture) ---
  p1 <- SpatialFeaturePlot(obj, features = target_cell_1, 
                           pt.size.factor = 2.4, 
                           alpha = c(0, 0.92), # Add a little bit of transparency to retain layers
                           stroke = 0) +
    
    # Key: values = c(0, 0.3, 1) makes gold appear earlier, making the image brighter
    scale_fill_gradientn(colours = c(col_trm_low, col_trm_mid, col_trm_high),
                         values = c(0, 0.3, 1)) + 
    
    ggtitle(paste(name, "CD8+ Trm")) +
    theme(legend.position = "right",
          plot.title = element_text(size=14, face="bold", color="black"))
  
  # --- 2. Plot Treg (Glacial texture) ---
  p2 <- SpatialFeaturePlot(obj, features = target_cell_2, 
                           pt.size.factor = 2.4, 
                           alpha = c(0, 0.92), 
                           stroke = 0) +
    
    scale_fill_gradientn(colours = c(col_treg_low, col_treg_mid, col_treg_high),
                         values = c(0, 0.3, 1)) +
    
    ggtitle(paste(name, "Th1 Program Treg")) +
    theme(legend.position = "right",
          plot.title = element_text(size=14, face="bold", color="black"))
  
  # Display combined plots
  print(p1 | p2)
}

#==============================================================================
message(">>> Starting Co-localization Score analysis...")

for (name in names(obj_list)) {
  obj <- obj_list[[name]]
  DefaultAssay(obj) <- "cell2location"
  
  # 1. Extract expression data
  # Note: Seurat reads usually go into the 'data' slot
  trm_expr <- GetAssayData(obj, slot = "data")[target_cell_1, ]
  treg_expr <- GetAssayData(obj, slot = "data")[target_cell_2, ]
  
  # 2. Calculate co-localization score (product)
  co_score <- trm_expr * treg_expr
  
  # 3. Add to Metadata for plotting
  obj <- AddMetaData(obj, metadata = co_score, col.name = "Trm_Treg_Co_Score")
  
  # Update the object in the list (This step is crucial, otherwise metadata won't be saved)
  obj_list[[name]] <- obj
  
  # 4. Plotting
  p <- SpatialFeaturePlot(obj, features = "Trm_Treg_Co_Score", 
                          pt.size.factor = 1.8,  # Slightly increase size
                          alpha = c(0.1, 1), 
                          stroke = 0) +
    scale_fill_viridis_c(option = "plasma") +    # Use plasma color scheme, highlighted areas are more prominent
    ggtitle(paste(name, ": Co-localization (Trm * Treg)"))
  
  print(p)
}

#===============================================================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales) # Used for normalizing data

# ================= 🔧 Manual RGB blending plotting function 🔧 =================
# This is a universal function, you can use it later to check colocalization of any two genes/cells
Plot_RGB_Overlay <- function(obj, feature_red, feature_green, sample_name) {
  
  # 1. Extract coordinates and expression data
  # Note: GetTissueCoordinates returns slightly different values in different Seurat versions
  # Compatible with Seurat v4/v5 here
  coords <- GetTissueCoordinates(obj, image = names(obj@images)[1])
  colnames(coords) <- c("imagerow", "imagecol") # Unify column names
  
  # Extract cell abundance data
  data <- FetchData(obj, vars = c(feature_red, feature_green))
  colnames(data) <- c("Red_Feature", "Green_Feature")
  
  # Combine
  plot_data <- cbind(coords, data)
  
  # 2. Key step: Data normalization (0-1)
  # cell2location scores vary greatly, must be compressed between 0-1 for color blending
  # Use quantile cutoff to prevent extreme values from making the whole image dark (similar to vmax='p99')
  limit_red <- quantile(plot_data$Red_Feature, 0.99)
  limit_green <- quantile(plot_data$Green_Feature, 0.99)
  
  plot_data$R_scaled <- pmin(plot_data$Red_Feature / limit_red, 1)
  plot_data$G_scaled <- pmin(plot_data$Green_Feature / limit_green, 1)
  
  # 3. Calculate blended color
  # Logic: R channel from Feature1, G channel from Feature2, B channel set to 0
  # alpha (transparency) determined by the max of both, transparent where neither is expressed
  plot_data$color_hex <- rgb(
    red = plot_data$R_scaled, 
    green = plot_data$G_scaled, 
    blue = 0, 
    alpha = pmax(plot_data$R_scaled, plot_data$G_scaled) * 0.9 # Transparency fine-tuning
  )
  
  # 4. Plotting
  # Note: ggplot's y-axis usually needs to be reversed to match the slice direction
  ggplot(plot_data, aes(x = imagecol, y = -imagerow)) +
    # Set background to black or dark color, simulating fluorescence microscope effect
    geom_point(color = plot_data$color_hex, size = 1.8, shape = 16) +
    theme_void() + # Remove axes
    theme(
      plot.background = element_rect(fill = "black"), # Black background
      panel.background = element_rect(fill = "black"),
      plot.title = element_text(color = "white", size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(color = "grey80", size = 10, hjust = 0.5)
    ) +
    labs(
      title = paste(sample_name, "Co-localization"),
      subtitle = paste0("Red: ", feature_red, " | Green: ", feature_green, " | Yellow: Co-localized")
    ) +
    # Force fixed aspect ratio to prevent slice distortion
    coord_fixed()
}

# ================= 🚀 Execute plotting 🚀 =================

message(">>> Starting to plot RGB colocalization (Red+Green=Yellow)...")

for (name in names(obj_list)) {
  obj <- obj_list[[name]]
  DefaultAssay(obj) <- "cell2location"
  
  tryCatch({
    p <- Plot_RGB_Overlay(obj, 
                          feature_red = target_cell_1,   # Trm -> Red
                          feature_green = target_cell_2, # Treg -> Green
                          sample_name = name)
    print(p)
  }, error = function(e) {
    message(paste("Error plotting sample", name, ":", e$message))
  })
}

#================================ TLS Display =======================================
library(Seurat)
library(ggplot2)
library(viridis)
# ================= 🔧 Single channel fluorescence plotting function (Dark version) 🔧 =================
# ================= 🔧 Deep Black Gold plotting function 🔧 =================
Plot_Deep_Black_Gold <- function(obj, feature_name, sample_name) {
  
  # 1. Extract coordinates
  if (packageVersion("Seurat") >= "5.0.0") {
    coords <- GetTissueCoordinates(obj, image = names(obj@images)[1])
  } else {
    coords <- GetTissueCoordinates(obj, image = names(obj@images)[1])
  }
  colnames(coords) <- c("imagerow", "imagecol")
  
  # 2. Extract data
  data <- FetchData(obj, vars = feature_name)
  colnames(data) <- "Score"
  
  # Combine
  plot_data <- cbind(coords, data)
  
  # 3. Plotting
  ggplot(plot_data, aes(x = imagecol, y = -imagerow, color = Score)) +
    geom_point(size = 1.5, shape = 16) + 
    
    # === 🎨 Key Modification: Extremely darken mid-values ===
    # Color logic:
    # 0% (Start) -> Black
    # 60% (Middle) -> Extremely deep dark gold (#332200), visually almost black with a slight gloss
    # 100% (End) -> Bright yellow
    # This way, only the highest scoring parts will truly "light up"
    scale_color_gradientn(
      colours = c("black", "#332200", "#FFFF00"), 
      values = c(0, 0.35, 1) 
    ) + 
    
    # === Dark theme ===
    theme_void() + 
    theme(
      plot.background = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      plot.title = element_text(color = "white", size = 15, face = "bold", hjust = 0.5),
      legend.background = element_rect(fill = "black"),
      legend.text = element_text(color = "white"),
      legend.title = element_text(color = "white")
    ) +
    labs(
      title = paste(sample_name, "TLS Score"),
      color = "Score"
    ) +
    coord_fixed()
}

# ================= 🚀 Execute plotting 🚀 =================

if (exists("test_data_T13")) {
  message(">>> Plotting darkened version of TLS score...")
  
  p <- Plot_Deep_Black_Gold(test_data_T13, 
                            feature_name = "TLS_Score1", 
                            sample_name = "T13")
  print(p)
  
} else {
  warning("Object test_data_T13 not found")
}

#============================== Correlation Analysis =======================================
library(Seurat)
library(ggplot2)
library(ggpubr) # Plotting package specifically used for adding statistical parameters

# ================= 1. Extract Data =================
# Assuming we are analyzing the T13 sample (you can change it to T5 or T14)
obj <- test_data_T14
DefaultAssay(obj) <- "cell2location"

# Define the two cell types you want to analyze
cell_x <- "Treg Th1 Program"  # X-axis
cell_y <- "CD8+ Trm"          # Y-axis

# Extract these two columns of data from the Seurat object
df_cor <- FetchData(obj, vars = c(cell_x, cell_y))

# Rename columns for easier coding
colnames(df_cor) <- c("Treg_Score", "Trm_Score")

# ================= 2. Data Filtering (Key Step) =================
# Count the number of spots before filtering
n_total <- nrow(df_cor)

# Strategy: Keep only spots where both cell scores are greater than 0 (Double positive Spots)
# This allows the regression line to reflect their quantitative relationship in the "co-existing region"
df_filtered <- subset(df_cor, Treg_Score > 0 & Trm_Score > 0)

n_keep <- nrow(df_filtered)
message(paste(">>> Data filtering complete:"))
message(paste("    Original number of Spots:", n_total))
message(paste("    Retained double positive Spots:", n_keep, "(", round(n_keep/n_total*100, 1), "%)"))
message("    Background spots with 0 expression for either have been filtered out.")

# ================= 3. Plot correlation scatter plot =================

# Use ggscatter to plot, a wrapper for ggplot2, perfect for publication-level statistical plots
p <- ggscatter(df_filtered, 
               x = "Treg_Score", 
               y = "Trm_Score", 
               
               # === Visual Adjustments ===
               color = "black",      # Point color
               fill = "lightgray",   # Point fill color (if shape 21)
               shape = 21,           # Use circles with borders
               size = 2.5,           # Point size
               alpha = 0.6,          # Transparency to prevent overlap
               
               # === Statistical Features ===
               add = "reg.line",                                         # Add regression line
               add.params = list(color = "red", fill = "pink", size = 1.5), # Set regression line to red
               conf.int = TRUE,                                          # Show confidence interval shadow
               
               # === Axis Labels ===
               xlab = paste(cell_x, "Score"), 
               ylab = paste(cell_y, "Score"),
               title = paste("Correlation in", "T13")
) +
  
  # === Add R value and P value ===
  stat_cor(method = "pearson",    # Use Pearson correlation (use "spearman" if skewed distribution)
           label.x.npc = "left",  # Label position: left
           label.y.npc = "top",   # Label position: top
           size = 5) +            # Font size
  
  # === Beautify Theme ===
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p)



# View results without filtering
library(Seurat)
library(ggplot2)
library(ggpubr) 

# ================= 1. Extract Data =================
# Assuming we are analyzing the T13 sample (you can change it to T5 or T14)
obj <- test_data_T5
DefaultAssay(obj) <- "cell2location"

# Define the two cell types you want to analyze
cell_x <- "Treg Th1 Program"  # X-axis
cell_y <- "CD8+ Trm"          # Y-axis

# Extract these two columns of data from the Seurat object
df_cor <- FetchData(obj, vars = c(cell_x, cell_y))

# Rename columns
colnames(df_cor) <- c("Treg_Score", "Trm_Score")

# ================= 2. Data Statistics (Modified: No filtering) =================
# Use all data directly
df_plot <- df_cor 

n_total <- nrow(df_plot)
message(paste(">>> Data preparation complete (no filtering):"))
message(paste("    Total Spots analyzed:", n_total))
message("    All spots retained, including background spots with 0 expression.")

# ================= 3. Plot correlation scatter plot =================

# Use ggscatter to plot
p <- ggscatter(df_plot, 
               x = "Treg_Score", 
               y = "Trm_Score", 
               
               # === Visual Adjustments ===
               color = "black",      
               fill = "lightgray",   
               shape = 21,           
               size = 2.0,           # Slightly smaller points because the number of points increased
               alpha = 0.5,          # Lower transparency to avoid excessive overlap at 0 values
               
               # === Statistical Features ===
               add = "reg.line",                                           
               add.params = list(color = "#FF4500", fill = "pink", size = 1), 
               conf.int = TRUE,                                           
               
               # === Axis Labels ===
               xlab = paste(cell_x, "Score"), 
               ylab = paste(cell_y, "Score"),
               title = paste("Correlation in T5 (All Spots)") # Title notes it includes all spots
) +
  
  # === Add R value and P value ===
  # Note: With a large number of 0 values, data may not fit a normal distribution,
  # If Pearson results seem too affected by 0 values, you can change method to "spearman"
  stat_cor(method = "pearson",    
           label.x.npc = "left",  
           label.y.npc = "top",   
           size = 5) +            
  
  # === Beautify Theme ===
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p)