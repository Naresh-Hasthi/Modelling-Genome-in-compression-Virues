
import matplotlib.pyplot as plt

def plot_compression(original_size, compressed_size):
    labels = ['Original', 'Compressed']
    sizes = [original_size, compressed_size]
    plt.bar(labels, sizes, color=['blue', 'green'])
    plt.ylabel('Size (bytes)')
    plt.title('Genome Compression Result')
    plt.savefig("output/compression_chart.png")
    plt.show()
