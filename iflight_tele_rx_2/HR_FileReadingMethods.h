#pragma once
#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <string>
#include <cstring>
#include <functional>

#if defined(_WIN32)
#include <Windows.h>
#elif defined(__linux__)
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#endif

#define NO_HEAD 0
#define NEW_ZY_HEAD 1
#define OLD_ZY_HEAD 2

const enum class filetype { ZenIns, SADAT };

#pragma pack(push,1)
//struct IQ_file_parameters {
//	double fc;
//	double fs;
//	double bw;
//	long long filesize = 0;
//	long long headsize = 0;
//};
// 参数结构体
struct IQ_file_parameters {
	double fc;          // 中心频率
	double fs;          // 采样率
	double bw;          // 带宽
	long long filesize = 0; // 文件大小
	long long headsize = 0; // 头大小
	int numChannels = 0;    // 声道数
	int bitsPerSample = 0;  // 样本位数
};


struct SADAT_filehead {
	char head[4];
	double fc;
	float bw;
	float fs;
	char time[20];
};


struct ZY_head_t {
	char head[4];
	uint32_t head_len;
	uint64_t fc;
	uint32_t bw;
	uint32_t fs;
	float level_corr;
};

struct ZY_old_head_t {
	int64_t fc;
	uint32_t fs;
	float level_corr;
};

union ZY_head_u {
	ZY_head_t new_head;
	ZY_old_head_t old_head;
};

struct wav_filehead {
	// RIFF头
	char chunkID[4];        // "RIFF"
	uint32_t chunkSize;     // 文件总大小-8
	char format[4];         // "WAVE"

	// fmt子块
	char subchunk1ID[4];    // "fmt "
	uint32_t subchunk1Size; // PCM为16
	uint16_t audioFormat;   // PCM为1
	uint16_t numChannels;   // 声道数
	uint32_t sampleRate;    // 采样率
	uint32_t byteRate;      // 每秒字节数
	uint16_t blockAlign;    // 每个样本的字节数
	uint16_t bitsPerSample; // 每个样本的位数

	// data子块
	char subchunk2ID[4];    // "data"
	uint32_t subchunk2Size; // 音频数据大小
};

// 完整的prof chunk结构体
struct IQS_ProfChunk {
	// chunk标识
	char chunk_id[4];            // "prof"
	uint32_t chunk_size;        // 格式块长度，小端

	// 数据起始偏移（块内偏移0）
	uint16_t file_number;       // 文件编号，大端
	uint8_t magic1;            // 0x8c
	uint8_t magic2;            // 0x22
	uint8_t magic3;            // 0x52
	uint8_t magic4;            // 0x9b
	uint8_t protocol_ver_major; // 0x00
	uint8_t protocol_ver_minor; // 0x02
	uint32_t api_version;      // API版本信息，大端

	// ... 可能有其他保留字段 ...
	uint8_t reserved[52];      // 从52到108字节的填充（根据文档偏移计算）

	// 从偏移108/块内64开始
	uint16_t struct_length;    // IQS_Profile + IQS_StreamInfo 结构体字节长度，大端

	// IQS_Profile 结构体字段
	double CenterFreq_Hz;
	double RefLevel_dBm;
	uint32_t DecimateFactor;
	int32_t RxPort;
	uint32_t BusTimeout_ms;
	int32_t TriggerSource;
	int32_t TriggerEdge;
	int32_t TriggerMode;
	uint64_t TriggerLength;
	int32_t TriggerOutMode;
	int32_t TriggerOutPulsePolarity;
	double TriggerLevel_dBm;
	double TriggerLevel_SafeTime;
	double TriggerDelay;
	double PreTriggerTime;
	int32_t TriggerTimerSync;
	double TriggerTimer_Period;
	uint8_t EnableReTrigger;
	double ReTrigger_Period;
	uint16_t ReTrigger_Count;
	int32_t DataFormat;
	int32_t GainStrategy;
	int32_t Preamplifier;
	uint8_t AnalogIFBWGrade;
	uint8_t IFGainGrade;
	uint8_t EnableDebugMode;
	int32_t ReferenceClockSource;
	double ReferenceClockFrequency;
	uint8_t EnableReferenceClockOut;
	int32_t SystemClockSource;
	double ExternalSystemClockFrequency;
	double NativeIQSampleRate_SPS;
	uint8_t EnableIFAGC;
	int8_t Atten;
	int32_t DCCancelerMode;
	int32_t QDCMode;
	float QDCIGain;
	float QDCQGain;
	float QDCPhaseComp;
	int8_t DCCIOffset;
	int8_t DCCQOffset;
	int32_t LOOptimization;
	double Bandwidth;

	// IQS_StreamInfo 结构体字段
	double IQSampleRate;
	uint64_t PacketCount;
	uint64_t StreamSamples;
	uint64_t StreamDataSize;
	uint32_t PacketSamples;
	uint32_t PacketDataSize;
	uint32_t GainParameter;

	// DeviceInfo 结构体
	uint16_t device_info_length;  // DeviceInfo结构体字节长度，大端
	uint64_t DeviceUID;
	uint16_t Model;
	uint16_t HardwareVersion;
	uint32_t MFWVersion;
	uint32_t FFWVersion;
	// ... DeviceInfo可能还有其他字段
};

#pragma pack(pop)




using std::endl;
using std::cout;
using std::vector;
using std::string;

typedef std::complex<float> T_COMPLEX_F;

inline string getFileExtension(const std::string& fileName) {
	size_t dotPos = fileName.find_last_of('.');
	if (dotPos != std::string::npos) {
		return fileName.substr(dotPos + 1);
	}
	return ""; // 如果没有找到扩展名，返回空字符串
}

inline std::string getFilenameWithoutExtension(const std::string& filepath) {
	// Find the last directory separator
	size_t lastSlash = filepath.find_last_of("/\\");

	// Extract filename with extension
	std::string filename = (lastSlash != std::string::npos) ?
		filepath.substr(lastSlash + 1) : filepath;

	// Remove extension
	size_t lastDot = filename.find_last_of('.');
	if (lastDot != std::string::npos) {
		return filename.substr(0, lastDot);
	}

	return filename;
}

inline void delete_existing_file(const char* filepath)
{
	std::ifstream fin(filepath, std::ios_base::out | std::ios_base::binary);
	if (fin.is_open()) {
		fin.close();
		remove(filepath);
	}
	return;
}

#if defined(_WIN32)
inline std::string printLastErrorMsg()
{
	LPSTR messageBuffer = nullptr;
	size_t size = FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
		NULL, GetLastError(), MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), (LPSTR)&messageBuffer, 0, NULL);
	std::string message(messageBuffer, size);
	return message;
}

inline vector<std::string> readDir(const std::string& dir)
{
	vector<std::string> res;

	std::string pattern(dir);
	pattern.append("\\*");
	WIN32_FIND_DATAA data;
	HANDLE hFind;
	//std::wstring tmp_w;
	if ((hFind = FindFirstFileA(pattern.c_str(), &data)) != INVALID_HANDLE_VALUE) {
		do {
			//tmp_w = data.cFileName;
			res.push_back(data.cFileName);
		} while (FindNextFileA(hFind, &data) != 0);
		FindClose(hFind);
	}

	return res;
}
#elif defined(__linux__)
inline std::vector<std::string> readDir(const std::string& path) {
	std::vector<std::string> filenames;
	DIR* dir;
	struct dirent* ent;

	if ((dir = opendir(path.c_str())) != nullptr) {
		while ((ent = readdir(dir)) != nullptr) {
			// 跳过 "." 和 ".." 目录
			if (strcmp(ent->d_name, ".") == 0 || strcmp(ent->d_name, "..") == 0) {
				continue;
			}

			// 获取文件名
			std::string filename = ent->d_name;
			filenames.push_back(filename);
		}
		closedir(dir);
	}
	else {
		perror("无法打开目录");
	}

	return filenames;
}
#endif

inline void parse_head(IQ_file_parameters& filehead, std::ifstream& fin, const char* filepath)
{
	filehead.fc = 0;
	filehead.fs = 0;
	filehead.bw = 0;
	filehead.filesize = 0;
	filehead.headsize = 0;
	filehead.numChannels = 0;
	filehead.bitsPerSample = 0;

	fin.seekg(0, std::ios_base::end);
	long long filesize = fin.tellg();
	fin.seekg(0, std::ios_base::beg);

	auto ext = getFileExtension(filepath);

	if (ext == "iq") {
		// 首先尝试以新版本读取
		char head_ch[4];
		fin.read(head_ch, 4);
		if (string(head_ch) == "iq") {
			fin.seekg(0, std::ios_base::beg);
			ZY_head_t tmp_head;
			fin.read((char*)&tmp_head, sizeof(ZY_head_t));
			/*printf("CenterFre: %.4f GHz, SampleRate : %.4f MHz, LevelCorrection : %.2f\n",
				tmp_head.fc / (float)1e9,
				tmp_head.fs / (float)1e6,
				tmp_head.level_corr);*/
			fin.seekg(tmp_head.head_len, std::ios_base::beg);   // 重设文件头位置
			filehead.fc = tmp_head.fc;
			filehead.fs = tmp_head.fs;
			filehead.bw = tmp_head.bw;
			filehead.filesize = filesize;
			filehead.headsize = tmp_head.head_len;
		}
		else {
			fin.seekg(0, std::ios_base::beg);
			ZY_old_head_t tmp_head;
			fin.read((char*)&tmp_head, sizeof(ZY_old_head_t));
			/*printf("CenterFre: %.4f GHz, SampleRate : %.4f MHz, LevelCorrection : %.2f\n",
				tmp_head.fc / (float)1e9,
				tmp_head.fs / (float)1e6,
				tmp_head.level_corr);*/
			fin.seekg(16, std::ios_base::beg);
			filehead.fc = tmp_head.fc;
			filehead.fs = tmp_head.fs;
			filehead.bw = tmp_head.fs / 1.28;
			filehead.filesize = filesize;
			filehead.headsize = 16;
		}
	}
	else if (ext == "sadat") {
		SADAT_filehead tmp_head;
		fin.read((char*)&tmp_head, 40/*sizeof(SADAT_filehead)*/);
		fin.seekg(40, std::ios_base::beg);
		filehead.fc = tmp_head.fc;
		filehead.fs = tmp_head.fs;
		filehead.bw = tmp_head.bw;
		filehead.filesize = filesize;
		filehead.headsize = 40;
	}
	else if (ext == "wav") {

		std::function<double(double)> tmp_byteswap_double =
			[](double x) {
			// 一个临时函数，用于Harogic .wav文件读取中需要的大小字节交换
			// 同时进行一字节后移，别问为什么要后移一字节，就是临时需要这么做
			uint64_t a;
			memcpy(&a, &x, 8);
			a >>= 8;		// 后移一字节
			uint8_t b[8];
			std::memcpy(b, &a, 8);
			std::reverse(b, b + 8);
			std::memcpy(&x, b, 8);
			return x;
			};

		char riff[4], wave[4];
		fin.read(riff, 4);
		uint32_t riffSize;
		fin.read(reinterpret_cast<char*>(&riffSize), 4);
		fin.read(wave, 4);

		if (std::strncmp(riff, "RIFF", 4) != 0 || std::strncmp(wave, "WAVE", 4) != 0) {
			std::cerr << "Error: Not a valid WAV file" << std::endl;
			return;
		}

		// 循环读取各个 chunk，直到找到 "fmt " 和 "data"
		bool fmt_found = false, data_found = false;
		uint32_t fmt_chunk_size = 0;
		uint32_t data_chunk_size = 0;
		uint64_t data_offset = 0;

		while (fin && !data_found) {
			char chunk_id[4];
			uint32_t chunk_size = 0;
			fin.read(chunk_id, 4);
			if (!fin) break;
			fin.read(reinterpret_cast<char*>(&chunk_size), 4);
			if (!fin) break;

			if (std::strncmp(chunk_id, "fmt ", 4) == 0) {
				fmt_found = true;
				std::vector<char> fmt_buf(chunk_size);
				fin.read(fmt_buf.data(), chunk_size);

				if (chunk_size >= 16) {
					filehead.numChannels = *reinterpret_cast<uint16_t*>(&fmt_buf[2]);
					filehead.fs = *reinterpret_cast<uint32_t*>(&fmt_buf[4]);
					filehead.bitsPerSample = *reinterpret_cast<uint16_t*>(&fmt_buf[14]);
				}
			}
			else if (std::strncmp(chunk_id, "data", 4) == 0) {
				data_found = true;
				data_chunk_size = chunk_size;
				data_offset = fin.tellg(); // data块数据起始
				break;
			}
			else if (std::strncmp(chunk_id, "prof", 4) == 0) {
				/*
				* Harogic wav文件读取较为复杂，需要用到msgpack...
				* 此处暂时读取出其中的中心频点，带宽暂时又采样率除以1.25计算
				*/
				IQS_ProfChunk tmpProf;
				int cur_loc = fin.tellg();
				fin.seekg(cur_loc - 4 - 4, std::ios_base::beg);
				fin.read((char*)&tmpProf, sizeof(IQS_ProfChunk));

				double cf = tmp_byteswap_double(tmpProf.CenterFreq_Hz);

				filehead.fc = cf;
				filehead.bw = filehead.fs / 1.25;

				fin.seekg(cur_loc + tmpProf.chunk_size, std::ios_base::beg);
			}
			else if (std::strncmp(chunk_id, "iq  ", 4) == 0) {
				int cur_loc = fin.tellg();
				uint32_t tmp_uint32;
				uint64_t tmp_uint64;
				float tmp_float;
				fin.read(reinterpret_cast<char*>(&tmp_uint64), 8); filehead.fc = tmp_uint64;
				fin.read(reinterpret_cast<char*>(&tmp_uint32), 4); filehead.bw = tmp_uint32;
				fin.read(reinterpret_cast<char*>(&tmp_uint32), 4); filehead.fs = tmp_uint32;
				fin.seekg(cur_loc + chunk_size, std::ios_base::beg);
			}
			else {
				// 跳过无关块
				fin.seekg(chunk_size, std::ios_base::cur);
			}

			// 按WAV标准，偶数对齐
			if (chunk_size % 2 == 1)
				fin.seekg(1, std::ios_base::cur);
		}

		if (!fmt_found || !data_found) {
			std::cerr << "Error: WAV missing fmt or data chunk\n";
			return;
		}

		filehead.headsize = data_offset;  // data块数据起始偏移
		filehead.filesize = filesize;

		// 把文件指针跳到数据区
		fin.seekg(filehead.headsize, std::ios_base::beg);
	}

	return;
}

inline void write_head(std::ofstream& fout, double fc, double fs, double bw, filetype type)
{
	if (type == filetype::ZenIns) {
		ZY_head_t tmp;
		memset(&tmp, 0, sizeof(ZY_head_t));
		memcpy(tmp.head, "iq", 2);
		tmp.fc = fc;
		tmp.fs = fs;
		tmp.bw = bw;
		tmp.head_len = 100;
		tmp.level_corr = 0.0;
		fout.write((char*)&tmp, sizeof(ZY_head_t));

		int byre_residual = 100 - sizeof(ZY_head_t);
		vector<uint8_t> zeros(byre_residual, 0);
		fout.write((char*)zeros.data(), byre_residual);

	}

	return;
}