#include "HR_FileReadingMethods.h"
#include <fftw3.h>

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <fstream>
#include <string>
#include <iomanip>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


std::vector<std::complex<double>> // 下变频后复数信号
commalgo_ddc // 基于最大保持法估计中心频率并下变频
(
	const std::vector<std::complex<double>>& data_in, // 待处理信号
	int fft_len // FFT长度
) 
{
    int total_len = data_in.size();
    int num_chunks = std::round(static_cast<double>(total_len) / fft_len);

    if (num_chunks < 2) {
        throw std::runtime_error("bad argument: insufficient data length for fft_len");
    }

    // Generate Hann window
    std::vector<double> win(fft_len);
    for (int i = 0; i < fft_len; ++i) {
        win[i] = 0.5 * (1.0 - std::cos(2.0 * M_PI * i / (fft_len - 1)));
    }

    // FFTW plan and buffer setup
    std::vector<fftw_complex> in_w(fft_len);
    std::vector<fftw_complex> out(fft_len);
    fftw_plan plan = fftw_plan_dft_1d(fft_len, in_w.data(), out.data(), FFTW_FORWARD, FFTW_ESTIMATE);

    // Helper lambda to compute FFT of a chunk with windowing
    auto compute_fft = [&](int start_idx, std::vector<std::complex<double>>& dest) {
        for (int i = 0; i < fft_len; ++i) {
            double w = win[i];
            auto val = data_in[start_idx + i];
            in_w[i][0] = w * val.real();
            in_w[i][1] = w * val.imag();
        }
        fftw_execute(plan);
        for (int i = 0; i < fft_len; ++i) {
            dest[i] = std::complex<double>(out[i][0], out[i][1]);
        }
        };

    // Process the first chunk
    std::vector<std::complex<double>> fft_maxhold(fft_len);
    std::vector<std::complex<double>> fft_tmp(fft_len);

    int num = num_chunks;
    int pos_b = (num - 2) * fft_len; // 0-based index

    compute_fft(pos_b, fft_maxhold);
    for (int i = 0; i < fft_len; ++i) {
        fft_maxhold[i] = std::abs(fft_maxhold[i]);
    }

    // Loop through remaining chunks
    num = num - 2;
    while (num > 0) {
        pos_b = (num - 1) * fft_len;
        compute_fft(pos_b, fft_tmp);
        for (int i = 0; i < fft_len; ++i) {
            double mag = std::abs(fft_tmp[i]);
            if (fft_maxhold[i].real() < mag) {
                fft_maxhold[i] = mag;
            }
        }
        num--;
    }

    fftw_destroy_plan(plan);

    // FFT shift and convert to dB scale (10 * log10)
    std::vector<double> fft_db(fft_len);
    int half = fft_len / 2;
    for (int i = 0; i < fft_len; ++i) {
        int src_idx = (i + half) % fft_len; // Equivalent to fftshift
        double val = std::abs(fft_maxhold[src_idx]);
        fft_db[i] = 10.0 * std::log10(val + 1e-12); // Prevent log(0)
    }

    // Find peak value and threshold (3dB bandwidth)
    double peak_val = fft_db[0];
    for (double v : fft_db) {
        if (v > peak_val) peak_val = v;
    }

    double snr = 3.0;//这里定义最大带宽保持门限
    double threshold = peak_val - snr;

    // Find indices above threshold
    std::vector<int> idx;
    for (int i = 0; i < fft_len; ++i) {
        if (fft_db[i] >= threshold) {
            idx.push_back(i);
        }
    }

    if (idx.empty()) {
        throw std::runtime_error("No frequencies found above the threshold.");
    }

    // Calculate center frequency f_mid based on shifted frequency vector F
    int idx_first = idx.front();
    int idx_last = idx.back();

    double f_beg = (static_cast<double>(idx_first) - fft_len / 2.0) / fft_len;
    double f_end = (static_cast<double>(idx_last) - fft_len / 2.0) / fft_len;
    double f_mid = (f_beg + f_end) / 2.0;

    // Digital Down Conversion (DDC)
    std::vector<std::complex<double>> data_out(total_len);
    for (int n = 0; n < total_len; ++n) {
        double angle = -2.0 * M_PI * f_mid * n;
        std::complex<double> exp_val(std::cos(angle), std::sin(angle));
        data_out[n] = data_in[n] * exp_val;
    }

    return data_out;
}


std::vector<std::complex<double>> // 滤波后的复数信号
commalgo_lpf // 低通滤波器 (复数信号)
(
	const std::vector<std::complex<double>>& data, // 输入信号数据
	double fs, // 采样率
    double f_cutoff, // 截止频率
    int filter_order = 64 //滤波器阶数 (默认 64，阶数越高频响越陡峭，但计算量越大)
)
{
    int M = filter_order;
    double fc_norm = f_cutoff / fs;

    // 设计 Hamming 窗 FIR 滤波器系数
    std::vector<double> h(M + 1, 0.0);
    double sum_h = 0.0;

    for (int n = 0; n <= M; ++n) {
        double m_val = n - M / 2.0;
        if (std::abs(m_val) < 1e-9) {
            h[n] = 2.0 * M_PI * fc_norm;
        }
        else {
            h[n] = std::sin(2.0 * M_PI * fc_norm * m_val) / m_val;
        }
        // Hamming 窗
        double w = 0.54 - 0.46 * std::cos(2.0 * M_PI * n / M);
        h[n] *= w;
        sum_h += h[n];
    }

    // 归一化系数（确保直流增益为 1）
    for (int n = 0; n <= M; ++n) {
        h[n] /= sum_h;
    }

    // 对复数数据进行卷积滤波
    size_t N = data.size();
    std::vector<std::complex<double>> filtered(N, std::complex<double>(0.0, 0.0));

    for (size_t i = 0; i < N; ++i) {
        std::complex<double> acc(0.0, 0.0);
        for (int j = 0; j <= M; ++j) {
            if (static_cast<int>(i) - j >= 0) {
                acc += data[i - j] * h[j];
            }
        }
        filtered[i] = acc;
    }

    return filtered;
}




std::vector<std::complex<double>> 
resampleLinear // 线性插值法对复数信号重采样
(
	const std::vector<std::complex<double>>& input, // 输入的复数序列
	double p, // 上采样因子
	double q // 下采样因子
)
{
    if (input.empty() || p <= 0 || q <= 0) {
        return {};
    }

    size_t in_size = input.size();
    // 计算重采样后的输出长度
    size_t out_size = static_cast<size_t>(std::round(static_cast<double>(in_size) * p / q));
    std::vector<std::complex<double>> output(out_size);

    for (size_t i = 0; i < out_size; ++i) {
        // 计算当前输出点在输入数组中的浮点索引位置
        double src_idx = static_cast<double>(i) * q / p;

        size_t idx = static_cast<size_t>(std::floor(src_idx));
        double frac = src_idx - idx;

        // 边界处理
        if (idx >= in_size - 1) {
            output[i] = input.back();
        }
        else {
            // 线性插值公式: v = v0 * (1 - frac) + v1 * frac
            std::complex<double> v0 = input[idx];
            std::complex<double> v1 = input[idx + 1];
            output[i] = v0 * (1.0 - frac) + v1 * frac;
        }
    }

    return output;
}


void 
burst_detector // 突发信号检测函数
(
	const std::vector<std::complex<double>>& data, // 输入复数信号数据
	int win_len, // 双滑窗单窗窗长
	double snr, // 突发高于底噪的比值阈值
	std::vector<int>& burst_pos, // 突发起点位置向量
	std::vector<int>& burst_len // 突发长度向量
) 
{
	int len = static_cast<int>(data.size());

	if (len < 2 * win_len) {
		throw std::runtime_error("bad input: data length is less than 2 * win_len");
	}

	// 信号瞬时能量 (data .* conj(data))
	std::vector<double> s_power(len);
	for (int i = 0; i < len; ++i) {
		s_power[i] = std::norm(data[i]); // 等效于 real^2 + imag^2
	}

	// 计算单窗功率（滑动窗口优化）
	int num_windows = len - win_len + 1;
	std::vector<double> window_power(num_windows, 0.0);

	double current_sum = 0.0;
	for (int i = 0; i < win_len; ++i) {
		current_sum += s_power[i];
	}
	window_power[0] = current_sum;

	for (int i = 1; i < num_windows; ++i) {
		current_sum = current_sum - s_power[i - 1] + s_power[i + win_len - 1];
		window_power[i] = current_sum;
	}

	for (int i = 0; i < num_windows; ++i) {
		window_power[i] /= win_len;
	}

	// 双滑窗比值（加入微小量 eps 防止分母为零）
	const double eps = 1e-16;
	int div_len = num_windows - win_len;
	std::vector<double> window_div(div_len);
	for (int i = 0; i < div_len; ++i) {
		window_div[i] = window_power[win_len + i] / (window_power[i] + eps);
	}

	// 清空输出
	burst_pos.clear();
	burst_len.clear();

	int burst_begin_pos = -1;

	// 寻找突发开始与结束
	for (int i = 0; i < div_len; ++i) {
		if (burst_begin_pos < 0 && window_div[i] > snr) {
			burst_begin_pos = i;
		}
		else if (burst_begin_pos >= 0 && window_div[i] < (1.0 / snr)) {
			int burst_end_pos = i + 2 * win_len;
			burst_pos.push_back(burst_begin_pos);
			burst_len.push_back(burst_end_pos - burst_begin_pos + 1);
			burst_begin_pos = -1;
		}
	}

	// 处理循环结束时仍未闭合的突发
	//if (burst_begin_pos >= 0) {
	//	int burst_end_pos = div_len + 2 * win_len;
	//	burst_pos.push_back(burst_begin_pos);
	//	//int actual_end = std::min(len, burst_end_pos);
	//	int actual_end = (std::min)(len, burst_end_pos);
	//	burst_len.push_back(actual_end - burst_begin_pos + 1);
	//}
}







namespace commalgo {
	using Complex = std::complex<double>;

	// 前向声明辅助函数
	std::vector<Complex> commalgo_lora_modulator(double fs, int SF, double BW, const std::string& way, int curr_send_sym);
	double commalgo_freq_par(const std::vector<Complex>& data);


	int // 返回估计的 SF 值，若未检测到有效 SF，则返回 -1
	commalgo_PAR_classifier // 基于峰均比的分类器
	(
		const std::vector<Complex>& data, // 待测的复数信号
		double fs // 采样率
	) 
	{
		int sf = -1;
		double BW = 812.5e3;

		int sf_list[] = { 5, 6, 7, 8 };
		int rep_counts[] = { 12, 14, 12, 12 };
		std::vector<double> par_vals(4, 0.0);

		for (int i = 0; i < 4; ++i) {
			int SF = sf_list[i];
			int rep = rep_counts[i];

			auto upchirp0 = commalgo_lora_modulator(fs, SF, BW, "up", 0);

			std::vector<Complex> upchirp;
			upchirp.reserve(upchirp0.size() * rep);
			for (int r = 0; r < rep; ++r) {
				upchirp.insert(upchirp.end(), upchirp0.begin(), upchirp0.end());
			}

			if (data.size() > upchirp.size()) {
				std::vector<Complex> tmp(upchirp.size());
				for (size_t j = 0; j < upchirp.size(); ++j) {
					tmp[j] = upchirp[j] * data[j];
				}
				par_vals[i] = commalgo_freq_par(tmp);
			}
			else {
				par_vals[i] = 0.0;
			}
		}

		auto max_it = std::max_element(par_vals.begin(), par_vals.end());
		double max_par = *max_it;
		int idx = std::distance(par_vals.begin(), max_it);

		double th = 100.0;
		if (max_par > th) {
			sf = idx + 5;
		}

		return sf;
	}


	
	double // 返回峰均比 (PAR) 值
	commalgo_freq_par // 计算频域峰均比 (PAR)
	(
		const std::vector<Complex>& data // 时域输入信号
	) 
	{
		int n = static_cast<int>(data.size());
		if (n == 0) return 0.0;

		// 使用 FFTW 分配内存
		fftw_complex* in = reinterpret_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * n));
		fftw_complex* out = reinterpret_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * n));

		// 复制输入数据
		for (int i = 0; i < n; ++i) {
			in[i][0] = data[i].real();
			in[i][1] = data[i].imag();
		}

		// 创建 FFTW 规划并执行正向变换
		fftw_plan plan = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(plan);

		// 计算幅值、峰值和均值
		// 注：MATLAB 代码中 fftshift 仅用于画图显示，峰值和均值在全频谱范围内平移不变，可直接计算。
		double sum_val = 0.0;
		double peak_val = 0.0;

		for (int i = 0; i < n; ++i) {
			double real = out[i][0];
			double imag = out[i][1];
			double mag = std::sqrt(real * real + imag * imag);

			sum_val += mag;
			if (mag > peak_val) {
				peak_val = mag;
			}
		}

		// 释放 FFTW 资源
		fftw_destroy_plan(plan);
		fftw_free(in);
		fftw_free(out);

		double mean_val = sum_val / static_cast<double>(n);
		const double eps = 1e-16;
		return peak_val / (mean_val + eps);
	}

	std::vector<Complex> // 返回生成的 LoRa 调制符号
	commalgo_lora_modulator // 生成一个 LoRa 调制符号
	(
		double fs, 
		int SF, 
		double BW, 
		const std::string& way, 
		int curr_send_sym
	) 
	{
		long long N = 1LL << SF;
		double Ts = static_cast<double>(N) / BW;
		double f0 = -BW / 2.0;
		double k = BW / Ts;

		std::vector<double> t;
		for (double ti = 0.0; ti < Ts - 1.0 / fs + 1e-9; ti += 1.0 / fs) {
			t.push_back(ti);
		}

		if (curr_send_sym < 0 || curr_send_sym >= N) {
			throw std::invalid_argument("非法参数：curr_send_sym 越界");
		}

		std::vector<Complex> lora(t.size());
		const double pi = acos(-1.0);

		if (way == "up") {
			for (size_t i = 0; i < t.size(); ++i) {
				double t_shifted = fmod(t[i] + curr_send_sym / BW, Ts);
				double phase = 2.0 * pi * (f0 * t_shifted + 0.5 * k * t_shifted * t_shifted);
				lora[i] = std::exp(Complex(0.0, phase));
			}
		}
		else if (way == "down") {
			for (size_t i = 0; i < t.size(); ++i) {
				double t_shifted = fmod(t[i] - curr_send_sym / BW, Ts);
				double phase = 2.0 * pi * (f0 * t_shifted + 0.5 * k * t_shifted * t_shifted);
				lora[i] = std::conj(std::exp(Complex(0.0, phase)));
			}
		}
		else {
			throw std::invalid_argument("非法参数：way 必须是 'up' 或 'down'");
		}

		return lora;
	}


	std::vector<double> // 返回归一化互相关函数值
	commalgo_cross_corr // 归一化互相关函数
	(
		const std::vector<Complex>& uw, // 独特码信号
		const std::vector<Complex>& data, // 待测信号
		int search_start, // 搜寻起点，输入1就行
		int search_end // 搜寻终点，不用太长，节约计算量
	) 
	{
		int m = static_cast<int>(uw.size());
		int n = static_cast<int>(data.size());

		if (search_end + m - 1 > n) {
			throw std::runtime_error("bad input: search range exceeds data length");
		}

		int cnt = search_end - search_start + 1;
		std::vector<double> corr(cnt, 0.0);

		double uw_energy = 0.0;
		for (int j = 0; j < m; ++j) {
			uw_energy += std::norm(uw[j]);
		}

		int search_data_len = search_end - search_start + m;
		int sd_start = search_start - 1;

		std::vector<double> search_data_energy(search_data_len);
		for (int i = 0; i < search_data_len; ++i) {
			search_data_energy[i] = std::norm(data[sd_start + i]) * uw_energy;
		}

		double R = 0.0;
		for (int i = 0; i < cnt; ++i) {
			int win_start = (search_start - 1) + i;
			Complex sum(0.0, 0.0);
			for (int j = 0; j < m; ++j) {
				sum += std::conj(uw[j]) * data[win_start + j];
			}
			double P_abs_sq = std::norm(sum); // |sum|^2

			if (i == 0) {
				R = 0.0;
				for (int j = 0; j < m; ++j) {
					R += search_data_energy[j];
				}
			}
			else {
				R = R - search_data_energy[i - 1] + search_data_energy[i + m - 1];
			}

			const double eps = 1e-16;
			corr[i] = P_abs_sq / (R + eps);
		}

		return corr;
	}


	int // 返回首峰索引，若未找到则返回 -1
	commalgo_find_first_peak // 寻首峰函数
	(
		const std::vector<double>& corr, // 归一化互相关函数图像
		double threshold = 0.7, // 阈值
		int start_idx = 2 // 起始索引2
	) 
	{
		int N = corr.size();
		int i = start_idx - 1; // 转换为 0-based 索引
		if (i < 1 || i >= N - 1 || N < 3) return -1;

		while (i < N - 1) {
			if ((corr[i] > threshold) && (corr[i] > corr[i - 1]) && (corr[i] >= corr[i + 1])) {
				return i;
			}
			i++;
		}
		return -1;
	}

	
	void 
	commalgo_cfo_correct // 小数倍 CFO 纠正
	(
		std::vector<Complex>& sig, // 待纠正信号
		int D, // 周期长度
		int idx // 起点索引
	) 
	{
		int len = sig.size();
		Complex P(0, 0);
		for (int i = 0; i < D; ++i) {
			P += std::conj(sig[idx + i]) * sig[idx + D + i];
		}
		double cfo_est = std::arg(P) / (2.0 * M_PI * D);
		for (int i = 0; i < len; ++i) {
			double angle = -2.0 * M_PI * cfo_est * i;
			sig[i] *= std::polar(1.0, angle);
		}
	}

	
	int // 返回解调符号值
	commalgo_lora_demodulator // 单个LoRa符号解调器
	(
		double fs, // 采样率
		int SF, // 扩频因子
		double BW, // 带宽
		int OSR, // 提高倍数
		const std::vector<Complex>& curr_sym // 当前符号数据
	) 
	{
		int N = 1 << SF;
		auto upchirp = commalgo_lora_modulator(fs, SF, BW, "up", 0);
		auto downchirp = commalgo_lora_modulator(fs, SF, BW, "down", 0);

		int sym_len = curr_sym.size();
		std::vector<Complex> dechirped1(sym_len);
		std::vector<Complex> dechirped2(sym_len);

		for (int i = 0; i < sym_len; ++i) {
			dechirped1[i] = curr_sym[i] * downchirp[i];
			dechirped2[i] = curr_sym[i] * upchirp[i];
		}

		// FFTW 规划与内存分配
		int fft_len = sym_len;
		fftw_complex* in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fft_len);
		fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fft_len);
		fftw_plan plan = fftw_plan_dft_1d(fft_len, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

		std::vector<Complex> fft_res1(fft_len), fft_res2(fft_len);

		// 处理 dechirped1 FFT
		for (int i = 0; i < fft_len; ++i) {
			in[i][0] = dechirped1[i].real();
			in[i][1] = dechirped1[i].imag();
		}
		fftw_execute(plan);
		for (int i = 0; i < fft_len; ++i) {
			fft_res1[i] = Complex(out[i][0], out[i][1]);
		}

		// 处理 dechirped2 FFT
		for (int i = 0; i < fft_len; ++i) {
			in[i][0] = dechirped2[i].real();
			in[i][1] = dechirped2[i].imag();
		}
		fftw_execute(plan);
		for (int i = 0; i < fft_len; ++i) {
			fft_res2[i] = Complex(out[i][0], out[i][1]);
		}

		fftw_destroy_plan(plan);
		fftw_free(in);
		fftw_free(out);

		// 累加 OSR 能量块
		std::vector<double> fft_mag1(N, 0.0);
		std::vector<double> fft_mag2(N, 0.0);

		for (int j = 0; j < OSR; ++j) {
			for (int i = 0; i < N; ++i) {
				fft_mag1[i] += std::abs(fft_res1[j * N + i]);
				fft_mag2[i] += std::abs(fft_res2[j * N + i]);
			}
		}

		auto max1_it = std::max_element(fft_mag1.begin(), fft_mag1.end());
		auto max2_it = std::max_element(fft_mag2.begin(), fft_mag2.end());

		int s1 = std::distance(fft_mag1.begin(), max1_it);
		int s2 = std::distance(fft_mag2.begin(), max2_it);

		return (*max1_it > *max2_it) ? s1 : s2;
	}

	
	bool // 返回是否成功解调 
	commalgo_lora_demod // 解调LoRa突发信号
	(
		std::vector<Complex>& data, // 输入数据
		double fs_osr, // 采样率
		int SF, // 扩频因子
		double BW, // 带宽
		int OSR, // 提高倍数
		int N_OSR, // 提高后的符号长度
		int corr_sym_num, // 互相关要用的符号数量
		int pilot_num, // 前导符号数量
		double pos_cfo, // 精同步字位置
		int demod_cnt, // 包含符号数量
		int pos_jump, // 0.25符号跳变位置
		std::vector<int>& demod_symbol // 输出解调符号
	) 
	{
		demod_symbol.clear();
		auto lora = commalgo_lora_modulator(fs_osr, SF, BW, "down", 0);

		// 第一次符号同步

		auto corr = commalgo_cross_corr(lora, data, 1, corr_sym_num * N_OSR);

		int first_peak_idx = commalgo_find_first_peak(corr, 0.5, 2);
		if (first_peak_idx == -1) return false;

		// 小数倍 CFO 纠正
		commalgo_cfo_correct(data, N_OSR, first_peak_idx);
		commalgo_cfo_correct(data, std::round(pilot_num / 4.0) * N_OSR, first_peak_idx);
		commalgo_cfo_correct(data, std::round(pilot_num / 2.0) * N_OSR, first_peak_idx);

		// 整数倍 CFO 纠正
		int sym_start = first_peak_idx + static_cast<int>(pos_cfo * N_OSR);

		if (sym_start < 0 ||
			sym_start > data.size() ||
			N_OSR > data.size() - sym_start) {
			return false;
		}

		std::vector<Complex> curr_sym(data.begin() + sym_start, data.begin() + sym_start + N_OSR);
		int symbol_u = commalgo_lora_demodulator(fs_osr, SF, BW, OSR, curr_sym);

		double myans = symbol_u;
		if (myans >= std::pow(2, SF) / 2.0) {
			myans -= std::pow(2, SF);
		}
		myans = myans / 2.0;

		double f_bin = BW / std::pow(2, SF);
		int data_len = data.size();
		for (int i = 0; i < data_len; ++i) {
			double angle = -2.0 * M_PI * myans * f_bin / fs_osr * i;
			data[i] *= std::polar(1.0, angle);
		}

		// 二次同步
		corr = commalgo_cross_corr(lora, data, 1, corr_sym_num * N_OSR);
		first_peak_idx = commalgo_find_first_peak(corr, 0.7, 2);
		if (first_peak_idx == -1) return false;

		int mv = 0;
		for (int i = 1; i <= demod_cnt; ++i) {
			if (i > pos_jump) {
				mv = N_OSR / 4;
			}
			int start = first_peak_idx + (i - 1) * N_OSR + mv;

			if (sym_start < 0 ||
				sym_start > data.size() ||
				N_OSR > data.size() - sym_start) {
				return false;
			}

			std::vector<Complex> c_sym(data.begin() + start, data.begin() + start + N_OSR);
			int curr_demod_symbol = commalgo_lora_demodulator(fs_osr, SF, BW, OSR, c_sym);
			demod_symbol.push_back(curr_demod_symbol);
		}

		return true;
	}


} // namespace commalgo








void main()
{
    // 尝试打开数据文件
    auto filepath_in = "C:\\Users\\11238\\Desktop\\iflight_tele_rx\\data\\lora\\mvlk500.wav";
    std::ifstream fin(filepath_in, std::ios_base::binary | std::ios_base::in);
    if (!fin.is_open()) {
        throw std::runtime_error("Unable to open input data file.");
    }
    // 读取接收数据
    IQ_file_parameters filehead;
    parse_head(filehead, fin, filepath_in);
    auto file_size = filehead.filesize;
    auto fs = filehead.fs;

    size_t num_samples = file_size / (2 * sizeof(int16_t));  // I+Q


    std::vector<std::complex<double>> iq(num_samples);

    int16_t i, q;
    for (size_t n = 0; n < num_samples; ++n) {
        fin.read(reinterpret_cast<char*>(&i), sizeof(int16_t));
        fin.read(reinterpret_cast<char*>(&q), sizeof(int16_t));
        iq[n] = std::complex<double>(
            i / 32768.0,
            q / 32768.0
        );
    }


#ifdef DEBUG
    auto data = iq;
    std::ofstream f("../log/seesee_iq.txt");
    int safe_len = 1000;
    for (size_t i = 0; i < min(size_t(safe_len), data.size()); ++i) {
        f << data[i].real() << " " << data[i].imag() << "\n";
    }
    std::cout << "dumped " << min(size_t(safe_len), data.size())
        << " samples to seesee_iq.txt\n";
#endif // DEBUG

    // 算法开始
    std::vector<int> burst_pos;
    std::vector<int> burst_len;
    burst_detector(iq, 1000, 10.0, /*win_len=1000, snr=10.0*/
        /*out*/ burst_pos, burst_len);

    // 循环处理
    for (size_t i = 0; i < burst_pos.size(); ++i) {
        std::cout << "Burst " << i + 1 << ": Position = " << burst_pos[i]
            << ", Length = " << burst_len[i] << "\n";
        // 临时拷贝（主要是我要对这段数据进行处理，直接在指针上面操作恐怕不是很好吧，因为会破坏原数据，但是原数据后续本身也会弃用，破坏了也不是不行）
        std::vector<std::complex<double>> burst(iq.begin() + burst_pos[i], iq.begin() + burst_pos[i] + burst_len[i]);
        // 基于最大保持法下变频
        std::vector<std::complex<double>>burst_ddc = commalgo_ddc(burst, 256);
        // 滤波
        std::vector<std::complex<double>> burst_ddc_fil = commalgo_lpf(burst_ddc, fs, 1e6);
        // 重采样
        double OSR = 16.0;
        double BW = 812.5e3;
        double fs_osr = BW * OSR;
        double p = fs_osr;
        double q = fs;
        std::vector<std::complex<double>> resampled_burst = resampleLinear(burst_ddc_fil, p, q);
        // debug 输出
        std::ofstream ofs("../log/resampled_burst.csv");
        ofs << "real,imag\n";
        for (const auto& v : resampled_burst) {
            ofs << v.real() << "," << v.imag() << "\n";
        }

        int sf = commalgo::commalgo_PAR_classifier(resampled_burst, fs_osr);

        // 解调
        int N_OSR = (1 << sf) * OSR;
        std::vector<int> demod_symbol; // 根据实际数据类型调整容器类型
        bool ok = false;

        N_OSR = (1 << sf) * OSR;
        demod_symbol.clear();
        double timelen = burst_len[i] / fs;

        if (sf == 5) {
            // 500
            if (timelen < 1.6e-3 && timelen > 1.4e-3) {
                ok = commalgo::commalgo_lora_demod(resampled_burst, fs_osr, sf, BW, OSR, N_OSR, 5, 12, 15, 38, 16, demod_symbol);
            }
            // 333
            if (timelen < 2.5e-3 && timelen > 2.3e-3) {
                ok = commalgo::commalgo_lora_demod(resampled_burst, fs_osr, sf, BW, OSR, N_OSR, 5, 12, 15, 60, 16, demod_symbol);
            }
        }
        else if (sf == 6) {
            // 250
            ok = commalgo::commalgo_lora_demod(resampled_burst, fs_osr, sf, BW, OSR, N_OSR, 10, 14, 17, 42, 18, demod_symbol);
        }
        else if (sf == 7) {
            // 150
            if (timelen < 6.0e-3 && timelen > 5.8e-3) {
                ok = commalgo::commalgo_lora_demod(resampled_burst, fs_osr, sf, BW, OSR, N_OSR, 5, 12, 15, 37, 16, demod_symbol);
            }
            // 100
            if (timelen < 7.7e-3 && timelen > 7.5e-3) {
                ok = commalgo::commalgo_lora_demod(resampled_burst, fs_osr, sf, BW, OSR, N_OSR, 5, 12, 15, 47/*48*/, 16, demod_symbol);
            }
        }
        else if (sf == 8) {
            // 50
            ok = commalgo::commalgo_lora_demod(resampled_burst, fs_osr, sf, BW, OSR, N_OSR, 5, 12, 15, 33/*这里改成33bug消失，之前是34.*/, 16, demod_symbol);
        }
        else {
            // sf = -1;
            // std::cout << "漏" << std::endl;
        }
		if (ok)
		{
			for (int val : demod_symbol) {
				std::cout << std::setw(4) << std::setfill(' ') << val;
			}
			std::cout << std::endl;
		}


    }

}