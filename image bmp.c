#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define WIDTH  256
#define HEIGHT 256

#pragma pack(push, 1) // Ensure structures are byte-aligned

typedef struct {
    uint16_t bfType;
    uint32_t bfSize;
    uint16_t bfReserved1;
    uint16_t bfReserved2;
    uint32_t bfOffBits;
} BMPFileHeader;

typedef struct {
    uint32_t biSize;
    int32_t  biWidth;
    int32_t  biHeight;
    uint16_t biPlanes;
    uint16_t biBitCount;
    uint32_t biCompression;
    uint32_t biSizeImage;
    int32_t  biXPelsPerMeter;
    int32_t  biYPelsPerMeter;
    uint32_t biClrUsed;
    uint32_t biClrImportant;
} BMPInfoHeader;

#pragma pack(pop)

void write_bmp(const char *filename, uint8_t image[HEIGHT][WIDTH][3])
{
    FILE *fp = fopen(filename, "wb");
    if (!fp)
    {
        perror("File opening failed");
        return;
    }

    int row_padded = (WIDTH * 3 + 3) & (~3); // BMP rows must be 4-byte aligned
    int data_size = row_padded * HEIGHT;
    int file_size = 54 + data_size;

    // BMP Header
    BMPFileHeader file_header = {0};
    file_header.bfType = 0x4D42; // 'BM'
    file_header.bfSize = file_size;
    file_header.bfOffBits = 54;

    // DIB Header
    BMPInfoHeader info_header = {0};
    info_header.biSize = 40;
    info_header.biWidth = WIDTH;
    info_header.biHeight = HEIGHT;
    info_header.biPlanes = 1;
    info_header.biBitCount = 24;
    info_header.biCompression = 0;
    info_header.biSizeImage = data_size;

    // Write headers
    fwrite(&file_header, sizeof(BMPFileHeader), 1, fp);
    fwrite(&info_header, sizeof(BMPInfoHeader), 1, fp);

    // Write pixel data (bottom-up order)
    uint8_t padding[3] = {0, 0, 0}; // Padding for alignment
    for (int y = HEIGHT - 1; y >= 0; y--)
    {
        fwrite(image[y], 3, WIDTH, fp);  // Write a row
        fwrite(padding, 1, row_padded - WIDTH * 3, fp); // Padding
    }

    fclose(fp);
    printf("BMP Image saved as %s\n", filename);
}

int main()
{
    uint8_t image[HEIGHT][WIDTH][3];

    // Generate a simple gradient pattern
    for (int y = 0; y < HEIGHT; y++)
    {
        for (int x = 0; x < WIDTH; x++)
        {
            image[y][x][0] = (uint8_t)x;   // Red
            image[y][x][1] = (uint8_t)y;   // Green
            image[y][x][2] = (uint8_t)(x + y) % 256; // Blue
        }
    }

    write_bmp("output.bmp", image);

    return 0;
}
