import streamlit as st
import matplotlib.pyplot as plt
import base64
from collections import Counter

from utils import (
    read_fasta, compute_entropy, compare_genomes, plot_codon_usage,
    gzip_compress, rle_compress, lzw_compress, bz2_compress, lzma_compress,
    huffman_compress, huffman_decompress
)

# -------------------------------
# Set custom background image + bold black text
# -------------------------------
def set_bg_image(image_path):
    with open(image_path, "rb") as img_file:
        encoded = base64.b64encode(img_file.read()).decode()
    st.markdown(
        f"""
        <style>
        .stApp {{
            background-image: url("data:image/jpg;base64,{encoded}");
            background-size: cover;
            background-attachment: fixed;
            background-repeat: no-repeat;
        }}

        html, body, [class*="st-"] {{
            color: black !important;
            font-weight: bold !important;
        }}

        .css-1d391kg, .css-1v0mbdj, .css-10trblm, .css-1cpxqw2 {{
            color: black !important;
            font-weight: bold !important;
        }}
        </style>
        """,
        unsafe_allow_html=True
    )

# ✅ Apply background (your actual image path)
set_bg_image(r"C:\Users\nares\OneDrive\Desktop\Biology\virus_genome_compression_project_fixed\Vir.jpg")

# -------------------------------
# Streamlit App Configuration
# -------------------------------
st.set_page_config(page_title="Virus Genome Compression Project", layout="wide")
st.title("🧬 Virus Genome Compression & Analysis Dashboard")

# Sidebar Uploads
st.sidebar.header("📂 Upload FASTA Files")
file1 = st.sidebar.file_uploader("Upload Genome 1", type=["fasta", "fa", "fna"])
file2 = st.sidebar.file_uploader("Upload Genome 2", type=["fasta", "fa", "fna"])

seq1 = seq2 = None
if file1 and file2:
    try:
        seq1 = read_fasta(file1)
        seq2 = read_fasta(file2)
    except Exception as e:
        st.sidebar.error(f"❌ Error reading uploaded files: {e}")

# Tabs
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "📊 Genome Stats", "🧬 Codon Usage", "📈 Entropy Plot",
    "📦 Compression Comparison", "🧠 Huffman Coding"
])

# -------------------------------
# TAB 1: Genome Stats
# -------------------------------
with tab1:
    st.header("📊 Genome Statistics")
    if seq1 and seq2:
        stats = compare_genomes(seq1, seq2)
        st.json(stats)
    else:
        st.warning("⚠️ Please upload both genomes.")

# -------------------------------
# TAB 2: Codon Usage
# -------------------------------
with tab2:
    st.header("🧬 Codon Usage Comparison")
    if seq1 and seq2:
        fig = plot_codon_usage(seq1, seq2)
        st.pyplot(fig)
    else:
        st.warning("⚠️ Please upload both genomes.")

# -------------------------------
# TAB 3: Entropy
# -------------------------------
with tab3:
    st.header("📈 Shannon Entropy (100 bp windows)")
    if seq1 and seq2:
        entropy1 = compute_entropy(seq1)
        entropy2 = compute_entropy(seq2)

        col1, col2 = st.columns(2)
        with col1:
            st.subheader("Genome 1")
            fig1, ax1 = plt.subplots()
            ax1.plot(entropy1, color="blue")
            ax1.set_xlabel("Window #")
            ax1.set_ylabel("Entropy")
            ax1.set_title("Genome 1 Entropy")
            st.pyplot(fig1)

        with col2:
            st.subheader("Genome 2")
            fig2, ax2 = plt.subplots()
            ax2.plot(entropy2, color="orange")
            ax2.set_xlabel("Window #")
            ax2.set_ylabel("Entropy")
            ax2.set_title("Genome 2 Entropy")
            st.pyplot(fig2)
    else:
        st.warning("⚠️ Upload both genomes to show entropy plots.")

# -------------------------------
# TAB 4: Compression
# -------------------------------
with tab4:
    st.header("📦 Compression Comparison")
    algo = st.selectbox("Choose Compression Algorithm", [
        "gzip", "Run-Length Encoding (RLE)", "LZW", "BZ2", "LZMA"
    ])

    if seq1 and seq2:
        size_chart_labels = []
        size_chart_values = []

        for i, seq in enumerate([seq1, seq2], start=1):
            st.subheader(f"Genome {i}")
            original_size = len(seq.encode())

            try:
                if algo == "gzip":
                    compressed = gzip_compress(seq)
                    compressed_size = len(compressed)
                elif algo == "Run-Length Encoding (RLE)":
                    compressed = rle_compress(seq)
                    compressed_size = len(compressed.encode())
                elif algo == "LZW":
                    compressed = lzw_compress(seq)
                    compressed_size = len(compressed) * 2
                elif algo == "BZ2":
                    compressed = bz2_compress(seq)
                    compressed_size = len(compressed)
                elif algo == "LZMA":
                    compressed = lzma_compress(seq)
                    compressed_size = len(compressed)
                else:
                    st.warning("❌ Unknown algorithm.")
                    continue

                ratio = compressed_size / original_size
                st.write(f"📄 Original Size: {original_size} bytes")
                st.write(f"📦 Compressed Size: {compressed_size} bytes")
                st.write(f"📉 Compression Ratio: {ratio:.2f}")

                size_chart_labels += [f"Original {i}", f"{algo} {i}"]
                size_chart_values += [original_size, compressed_size]

            except Exception as e:
                st.error(f"❌ Error in Genome {i}: {e}")

        if len(size_chart_values) == 4:
            fig, ax = plt.subplots()
            ax.bar(size_chart_labels, size_chart_values, color=["gray", "green", "gray", "blue"])
            ax.set_ylabel("Size (bytes)")
            st.pyplot(fig)
    else:
        st.warning("⚠️ Please upload both genome files.")

# -------------------------------
# TAB 5: Huffman Compression
# -------------------------------
with tab5:
    st.header("🧠 Huffman Compression")

    if seq1 and seq2:
        try:
            enc1, codes1, bits1, compressed1 = huffman_compress(seq1)
            enc2, codes2, bits2, compressed2 = huffman_compress(seq2)

            st.subheader("Genome 1")
            st.write(f"Original: {bits1} bits, Compressed: {compressed1} bits, Ratio: {compressed1 / bits1:.2f}")

            st.subheader("Genome 2")
            st.write(f"Original: {bits2} bits, Compressed: {compressed2} bits, Ratio: {compressed2 / bits2:.2f}")

            if st.checkbox("🔍 Show Huffman Code Maps"):
                st.text("Genome 1 Code Map")
                st.json(codes1)
                st.text("Genome 2 Code Map")
                st.json(codes2)

            if st.checkbox("🔍 Show Encoded Sequences (first 500 bits)"):
                st.text("Genome 1")
                st.code(enc1[:500] + '...')
                st.text("Genome 2")
                st.code(enc2[:500] + '...')

            # Frequency Plot
            st.subheader("📊 Nucleotide Frequencies")
            freq1 = Counter(seq1)
            freq2 = Counter(seq2)
            nucs = sorted(set(freq1.keys()) | set(freq2.keys()))
            val1 = [freq1.get(n, 0) for n in nucs]
            val2 = [freq2.get(n, 0) for n in nucs]

            fig1, ax1 = plt.subplots()
            width = 0.35
            ax1.bar(range(len(nucs)), val1, width=width, label="Genome 1", color="blue")
            ax1.bar([x + width for x in range(len(nucs))], val2, width=width, label="Genome 2", color="orange")
            ax1.set_xticks([x + width / 2 for x in range(len(nucs))])
            ax1.set_xticklabels(nucs)
            ax1.set_ylabel("Count")
            ax1.set_title("Nucleotide Frequencies")
            ax1.legend()
            st.pyplot(fig1)

            # Code Length Plot
            st.subheader("📏 Huffman Code Lengths")
            len1 = [len(codes1.get(n, '')) for n in nucs]
            len2 = [len(codes2.get(n, '')) for n in nucs]

            fig2, ax2 = plt.subplots()
            ax2.bar(range(len(nucs)), len1, width=width, label="Genome 1", color="green")
            ax2.bar([x + width for x in range(len(nucs))], len2, width=width, label="Genome 2", color="red")
            ax2.set_xticks([x + width / 2 for x in range(len(nucs))])
            ax2.set_xticklabels(nucs)
            ax2.set_ylabel("Bits")
            ax2.set_title("Huffman Code Lengths")
            ax2.legend()
            st.pyplot(fig2)

        except Exception as e:
            st.error(f"❌ Huffman Error: {e}")
    else:
        st.warning("⚠️ Please upload both genome files.")

