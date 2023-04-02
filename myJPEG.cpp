#include<iostream>
#include"myJPEG.h"
#include"qdbmp.h"
#include<map>
#include<utility>
#define _USE_MATH_DEFINES
#include<math.h>


using namespace std;

struct {
    int height;
    int width;
} image;

struct {
    unsigned char id;
    unsigned char width;
    unsigned char height;
    unsigned char quant;
} subVector[4];

unsigned char maxWidth, maxHeight;

struct acCode {
    unsigned char len;
    unsigned char zeros;
    int value;
};

struct RGB {
    unsigned char R, G, B;
};

typedef double BLOCK[8][8];

const int DC = 0;
const int AC = 1;

unsigned int color_channel = 3;

int color_huff[4][2];

//Declare map of huffman tree (pair) and source code (unsigned char)
std::map<std::pair<unsigned char, unsigned int>, unsigned char> huffTable[2][2];

int quantTable[4][128];

double cos_cache[200];

void init_cos_cache() {
    for (int i = 0; i < 200; i++) {
        cos_cache[i] = cos(i * M_PI / 16.0);
    }
}

class MCU {
public:
    BLOCK mcu[4][2][2];
    void show() {
        showSectionHeader("MCU SHOW");
        for (int id = 1; id <= color_channel; id++) {
            for (int h = 0; h < subVector[id].height; h++) {
                for (int w = 0; w < subVector[id].width; w++) {
                    printf("mcu id: %d, %d %d\n", id, h, w);
                    for (int i = 0; i < 8; i++) {
                        for (int j = 0; j < 8; j++) {
                            printf("%lf ", mcu[id][h][w][i][j]);
                        }
                        printf("\n");
                    }
                }
            }
        }
    };
    void quantify() {
        for (int id = 1; id <= color_channel; id++) {
            for (int h = 0; h < subVector[id].height; h++) {
                for (int w = 0; w < subVector[id].width; w++) {
                    for (int i = 0; i < 8; i++) {
                        for (int j = 0; j < 8; j++) {
                            mcu[id][h][w][i][j] *= quantTable[subVector[id].quant][i * 8 + j];
                        }
                    }
                }
            }
        }
    };
    void zigzag() {
        for (int id = 1; id <= 3; id++) {
            for (int h = 0; h < subVector[id].height; h++) {
                for (int w = 0; w < subVector[id].width; w++) {
                    int zz[8][8] = {
                            { 0,  1,  5,  6, 14, 15, 27, 28},
                            { 2,  4,  7, 13, 16, 26, 29, 42},
                            { 3,  8, 12, 17, 25, 30, 41, 43},
                            { 9, 11, 18, 24, 31, 40, 44, 53},
                            {10, 19, 23, 32, 39, 45, 52, 54},
                            {20, 22, 33, 38, 46, 51, 55, 60},
                            {21, 34, 37, 47, 50, 56, 59, 61},
                            {35, 36, 48, 49, 57, 58, 62, 63}
                    };
                    for (int i = 0; i < 8; i++) {
                        for (int j = 0; j < 8; j++) {
                            zz[i][j] = mcu[id][h][w][zz[i][j] / 8][zz[i][j] % 8];
                        }
                    }
                    for (int i = 0; i < 8; i++) {
                        for (int j = 0; j < 8; j++) {
                            mcu[id][h][w][i][j] = zz[i][j];
                        }
                    }
                }
            }
        }
    };
    void idct() {
        for (int id = 1; id <= 3; id++) {
            for (int h = 0; h < subVector[id].height; h++) {
                for (int w = 0; w < subVector[id].width; w++) {
                    double tmp[8][8] = { 0 };
                    double s[8][8] = {};
                    for (int j = 0; j < 8; j++) {
                        for (int x = 0; x < 8; x++) {
                            for (int y = 0; y < 8; y++) {
                                s[j][x] += c(y) * mcu[id][h][w][x][y] * cos_cache[(j + j + 1) * y];
                            }
                            s[j][x] = s[j][x] / 2.0;
                        }
                    }
                    for (int i = 0; i < 8; i++) {
                        for (int j = 0; j < 8; j++) {
                            for (int x = 0; x < 8; x++) {
                                tmp[i][j] += c(x) * s[j][x] * cos_cache[(i + i + 1) * x];
                            }
                            tmp[i][j] = tmp[i][j] / 2.0;
                        }
                    }
                    for (int i = 0; i < 8; i++) {
                        for (int j = 0; j < 8; j++) {
                            mcu[id][h][w][i][j] = tmp[i][j];
                        }
                    }
                }
            }
        }
    }
    void decode() {
        this->quantify();
        //cout << "quantify\n";
        //this->show();
        this->zigzag();
        //cout << "zigzag\n";
        //this->show();
        this->idct();
        //cout << "idct\n";
        //this->show();
    }
    RGB** toRGB() {
        RGB** ret = (RGB**)malloc(sizeof(RGB**) * maxHeight * 8);
        for (int i = 0; i < maxHeight * 8; i++) {
            ret[i] = (RGB*)malloc(sizeof(RGB*) * maxWidth * 8);
        }
        double Y=0, Cb=0, Cr = 0;
        for (int i = 0; i < maxHeight * 8; i++) {
            for (int j = 0; j < maxWidth * 8; j++) {
                Y = trans(1, i, j);
                //printf("Y:%d \n", Y);
                Cb = trans(2, i, j);
                //printf("Cb:%d \n", Cb);
                Cr = trans(3, i, j);
                //printf("Cr:%d \n", Cr);

                ret[i][j].R = chomp(Y + 1.402 * Cr + 128);
                ret[i][j].G = chomp(Y - 0.34414 * Cb - 0.71414 * Cr + 128);
                ret[i][j].B = chomp(Y + 1.772 * Cb + 128);
            }
        }
        return ret;
    }

    private:
        double cc(int i, int j) {
            if (i == 0 && j == 0) {
                return 1.0 / 2.0;
            }
            else if (i == 0 || j == 0) {
                return 1.0 / sqrt(2.0);
            }
            else {
                return 1.0;
            }
        }
        double c(int i) {
            static double x = 1.0 / sqrt(2.0);
            if (i == 0) {
                return x;
            }
            else {
                return 1.0;
            }
        }
        unsigned char chomp(double x) {
            if (x > 255.0) {
                return 255;
            }
            else if (x < 0) {
                return 0;
            }
            else {
                return (unsigned char)x;
            }
        }
        double trans(int id, int h, int w) {          
            int vh = h * subVector[id].height / maxHeight;
            int vw = w * subVector[id].width / maxWidth;
            //printf("%d ", id);
            //printf("%d ", vh);
            //printf("%d ", vw);
            //printf("%d \n", mcu[id][vh / 8][vw / 8][vh % 8][vw % 8]);
            return mcu[id][vh / 8][vw / 8][vh % 8][vw % 8];
        }
};

//create huffman code tree. Pair len with code
//retun a pointer to pair 
std::pair<unsigned char, unsigned int>* createHuffCode(unsigned char* a, unsigned int number) 
{
    int si = sizeof(std::pair<unsigned char, unsigned int>);//算容器的内存
    auto ret = (std::pair<unsigned char, unsigned int>*)malloc(si * number);//智能指针 分配HuffTable的内存
    int code = 0;
    int count = 0;
    //
    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < a[i]; j++) {
            ret[count++] = std::make_pair(i + 1, code);//i+1是huffman tree叶子节点高度，code是编码值
            code += 1;
        }
        code = code << 1;
    }
    return ret;
}

void showSectionHeader(const char *s)
{
    cout<<"******"<<s<<"******\n";
}

unsigned int getSectionLength(unsigned char* &buf)
{
    unsigned int length;
    buf++;
    length = (*buf);
    buf++;
    length = length * 256 + (*buf);
    cout << "Section length: " << length << endl;
    return length;
}

void ParseDQT(unsigned char*& buf)
{
    showSectionHeader("DQT");
    unsigned int length;
    length = getSectionLength(buf);
    buf++;
    unsigned precision = (*buf) >> 4 == 0 ? 8 : 16;
    cout << "Precision: " << precision << endl;
    precision /= 8;
    unsigned char id = *(buf) & 0x0F;
    cout << "DQT ID: " << unsigned int(id) << endl;

    for (int i = 0; i < 64; i++) {
        unsigned char t = 0;
        for (int p = 0; p < precision; p++) {
            unsigned char s;
            buf++;
            s = *buf;
            t == t << 8;
            t += s;
        }
        quantTable[id][i] = t;
    }

    for (int i = 0; i < 64; i++) {
        if (i % 8 == 0) {
            printf("\n");
        }
        //printf("%2d ", quantTable[id][i]);
        cout << quantTable[id][i] << ' ';
    }
    cout << endl;

}

void ParseDHT(unsigned char *&buf)
{
    showSectionHeader("DHT");
    unsigned int length;
    length=getSectionLength(buf);
    buf++;
    unsigned char v[1];
    v[0] = *buf;
    unsigned char DCorAC=v[0] >> 4;
    unsigned char id = v[0] & 0x0F;
    cout<<"DHT ID: " << (DCorAC == 0 ? "DC " : "AC ")  <<unsigned int(id) << endl;

    unsigned char a[16];
    unsigned int number = 0;
    for (int i = 0; i < 16; i++) {
        buf++;
        a[i] = *buf;
        cout << unsigned int(a[i]) << " ";
        number += *buf;
    }
    cout << "Coded_Nums： " << number << endl;
    
    auto huffCode = createHuffCode(a, number);
    for (int i = 0; i < number; i++) {
        buf++;
        unsigned char v;
        v = *buf;
        //map of huffman tree and source code
        huffTable[DCorAC][id][huffCode[i]] = v;
        printf("%d  %d : %d\n", huffCode[i].first, huffCode[i].second, v);
    }
    free(huffCode);
}

void ParseSOF0(unsigned char *buf)
{
    showSectionHeader("SOF0");
    unsigned int length;
    length = getSectionLength(buf);
    buf++;
    unsigned precision = (*buf) >> 4 == 0 ? 8 : 16;
    cout << "Precision: " << precision << endl;
    precision /= 8;
    buf++;
    image.height = *(buf) * 256 + *(buf + 1);
    buf += 2;
    image.width = *(buf) * 256 + *(buf + 1);
    buf += 2;
    cout << "Image Size: " << image.height <<"x"<< image.width << endl;
    color_channel = *(buf);
    buf += 1;
    for (int i = 0; i < color_channel;i++)
    {
        cout << "Color id: " << int(buf[0]) << endl;
        cout << "Horizontal Sample " << int(buf[1]>>4) << endl;
        cout << "Vertical Sample: " << int(buf[1]&(0x0F)) << endl;
        cout << "QuantTable id: " << int(buf[2]) << endl;
        subVector[buf[0]].id = buf[0];
        subVector[buf[0]].width = buf[1]>>4;
        subVector[buf[0]].height = buf[1]&0x0F;
        subVector[buf[0]].quant = buf[2];
        maxHeight = (maxHeight > subVector[buf[0]].height ? maxHeight : subVector[buf[0]].height);
        maxWidth = (maxWidth > subVector[buf[0]].width ? maxWidth : subVector[buf[0]].width);
        buf += 3;
    } 
    buf += 3;
}

void ParseSOS(unsigned char*& buf)
{
    showSectionHeader("SOS");
    unsigned int length;
    length = getSectionLength(buf);
    buf++;
    color_channel = (*buf);
    buf++;
    for (int i = 0; i < color_channel; i++)
    {
        cout << "Color id: " << int(buf[0]) << endl;
        cout <<"DC_id " <<(int(buf[1] >> 4)) 
            <<"\nAC_id " << int(buf[1] & 0x0F) << endl;
        color_huff[int(buf[0])][DC] = buf[1] >> 4;
        color_huff[int(buf[0])][AC] = buf[1] &0x0F;

        buf += 2;
    }
    buf += 3;
}

bool getBit(unsigned char*&buf)
{
    static unsigned char buffer;
    static unsigned char count = 0;
    if (count == 0)
    {
        //buf--;
        buffer = *buf;
        buf += 1;
        //cout <<(int) buffer << endl;
        //printf("read %X\n", buffer);
        if (buffer == 0xFF)
        {
            unsigned char check;
            check = *(buf);
            buf += 1;
            if (check != 0x00)
            {
                fprintf(stderr, "data段错误");
            }
        }
    }
    //"&" operation to get bit
    bool ret = buffer & (1 << (7 - count));
    //count to go the next byte;
    //cout << int(ret) << ' ';
    count = (count == 7 ? 0 : count + 1);
    //cout << (count == 0 ? "\n" : "");
    
    
    return ret;
}

unsigned char matchHuff(unsigned char*&buf, unsigned char number, unsigned char ACorDC) {
    unsigned int len = 0;
    unsigned char codeLen;
    for (int count = 1; ; count++) {
        len = len << 1;
        //buf--;
        len += (unsigned int)getBit(buf);
        //cout << len << endl;
        //find huffman code (pair) from huffTable (map)
        if (huffTable[ACorDC][number].find(std::make_pair(count, len)) != huffTable[ACorDC][number].end()) {
            codeLen = huffTable[ACorDC][number][std::make_pair(count, len)];
            //codeLen is the source code
            //cout << int(codeLen) << endl;
            return codeLen;
        }
        if (count > 16) 
        { 
            cout<<"Do not find match key\n"; 
            count = 1; 
            len = 0; 
        }
    }
}

int readDC(unsigned char * &buf, unsigned char number) {
    unsigned char codeLen = matchHuff(buf, number, DC);
    //cout << int(codeLen) << endl;
    if (codeLen == 0) { return 0; }
    unsigned char first = getBit(buf);
    int ret = 1;
    // DC table decode
    for (int i = 1; i < codeLen; i++) {
        unsigned char b = getBit(buf);
        ret = ret << 1;
        ret += first ? b : !b;
    }
    ret = first ? ret : -ret;
    //printf("read DC: len %d, value %d\n", codeLen, ret);
    return ret;
}

// 計算ZRL
acCode readAC(unsigned char *&buf, unsigned char number) {
    unsigned char x = matchHuff(buf, number, AC);
    unsigned char zeros = x >> 4;
    unsigned char codeLen = x & 0x0F;
    if (x == 0) {
        return acCode{ 0,0,0 };
    }
    else if (x == 0xF0) {
        return acCode{ 0, 16, 0 };
    }
    unsigned  char first = getBit(buf);
    int code = 1;
    for (int i = 1; i < codeLen; i++) {
        unsigned char b = getBit(buf);
        code = code << 1;
        code += first ? b : !b;
    }
    code = first ? code : -code;
    //printf("read AC: %d %d %d\n", codeLen, zeros, code);
    return acCode{ codeLen, zeros, code };
}

MCU readMCU(unsigned char*& buf)
{
    static int dc[4] = { 0,0,0,0 };
    auto mcu = MCU();
    for (int i = 1; i <= color_channel; i++)
    {
        for (int h = 0; h < subVector[i].height; h++)
        {
            for (int w = 0; w < subVector[i].width; w++)
            {
                dc[i] = readDC(buf, i/2) + dc[i];
                mcu.mcu[i][h][w][0][0] = dc[i];
                unsigned int count = 1;
                while (count < 64) {
                    acCode ac = readAC(buf, i/2);
                    if (ac.len == 0 && ac.zeros == 16) {
                        for (int j = 0; j < ac.zeros; j++) {
                            mcu.mcu[i][h][w][count / 8][count % 8] = 0;
                            count++;
                        }
                    }
                    else if (ac.len == 0) {
                        break;
                    }
                    else {
                        for (int j = 0; j < ac.zeros; j++) {
                            mcu.mcu[i][h][w][count / 8][count % 8] = 0;
                            count++;
                        }
                        mcu.mcu[i][h][w][count / 8][count % 8] = ac.value;
                        count++;
                    }
                }
                while (count < 64) {
                    mcu.mcu[i][h][w][count / 8][count % 8] = 0;
                    count++;
                }
            }
        }
    }
    return mcu;
}

void ParseData(unsigned char* &buf, const char* outfilename)
{
    showSectionHeader("ReadData");
    int w = (image.width - 1) / (8 * maxWidth) + 1;
    int h = (image.height - 1) / (8 * maxHeight) + 1;

    BMP* bmp = BMP_Create(maxWidth * 8 * w, maxHeight * 8 * h, 24);
    
    for (int i = 0; i < h; i++)
    {
        for (int j = 0; j < w; j++)
        {
            MCU mcu = readMCU(buf);
            //mcu.show();
            mcu.decode();
            //mcu.show();
            RGB** b = mcu.toRGB();
            for (int y = i * 8 * maxHeight; y < (i + 1) * 8 * maxHeight; y++) {
                for (int x = j * 8 * maxWidth; x < (j + 1) * 8 * maxWidth; x++) {
                    int by = y - i * 8 * maxHeight;
                    int bx = x - j * 8 * maxWidth;
                    BMP_SetPixelRGB(bmp, x, y, b[by][bx].R, b[by][bx].G, b[by][bx].B);
                }
            }
        }
    }
    cout<<"Ready to output file\n";
    BMP_WriteFile(bmp, outfilename);
}

void Jpeg_Hexdecode(const char* infilename, const char* outfilename)
{
    
    init_cos_cache();
    FILE* fp;
    unsigned int length_of_file;
    //unsigned int width;
    //unsigned int height;
    unsigned char* buf;
    //unsigned char* components[3];

    /* Load the Jpeg into memory 将JPEG文件读入缓冲区*/
    fp = fopen(infilename, "rb");
    cout <<"Load filename : "<< infilename << endl;
    if (fp == NULL)
    {
        cout<<"Cannot open filename\n";
        return;
    }
    fseek(fp, 0, SEEK_END);
    length_of_file = ftell(fp);
    fseek(fp, 0, SEEK_SET);

    buf = (unsigned char*)malloc(length_of_file+4);
    unsigned char * buf_pointer=new unsigned char;
    buf_pointer = buf;
    if (buf == NULL)
    {
        cout<<"Not enough memory for loading file\n";
        return;
    }
    fread(buf, length_of_file, 1, fp);
    fclose(fp);
    cout << "File length : " << length_of_file << "bytes\n";
    for(int i=0;i<=length_of_file;i++)
    {
        if (*buf_pointer == 0xFF)
        {
            buf_pointer = buf_pointer + 1;
            switch (*buf_pointer)
            {
            case SOI_MARKER:
                //cout << "Start of Image\n";
                break;
            case APP0_MARKER:
            case 0xE1:
            case 0xE2:
            case 0xE3:
            case 0xE4:
            case 0xE5:
            case 0xE6:
            case 0xE7:
            case 0xE8:
            case 0xE9:
            case 0xEA:
            case 0xEB:
            case 0xEC:
            case 0xED:
            case 0xEE:
            case 0xEF:
                //cout << "APP segement\n";
                //readAPP(f);
                break;
            case COM_MARKER:
                //cout << "Comment\n";
                //readCOM(f);
                break;
            case DQT_MARKER:
                //cout << "Define Quantization Tables\n";
                ParseDQT(buf_pointer);
                break;
            case SOF0_MARKER:
                //cout << "Start Of Frame (Baseline)\n";
                ParseSOF0(buf_pointer);
                break;
            case SOF2_MARKER:
                //cout << "Start Of Frame (Progressive)\n";
                //readSOF2(f);
                break;
            case DHT_MARKER:
                //cout << "Define Huffman Table\n";
                ParseDHT(buf_pointer);
                break;
            case SOS_MARKER:
                //cout << "Start Of Scan\n";
                ParseSOS(buf_pointer);
                ParseData(buf_pointer, outfilename);
                goto label;
                //readData(f);
                break;
            case EOI_MARKER:
                cout << "End of Image\n";
                break;
            }
        }
        buf_pointer++;
    }
    label:
    free(buf);
    buf_pointer = NULL;
    delete buf_pointer;
}