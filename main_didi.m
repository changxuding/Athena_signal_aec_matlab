clear all; 
clc;

a = 'dt';
b = 'zhoujielun';
c = '/out.wav';
echo_file = ['./', a, '/', b, '/echo.wav'];
far_file = ['./', a, '/', b, '/far.wav'];
out_file = ['./', a, '/', b, c]; 

[y_echo, fs_echo] = audioread(echo_file,'native');
[y_far, fs_far] = audioread(far_file,'native');
y_echo = double(y_echo);
y_far = double(y_far);
y_out = zeros(1, length(y_echo));

% static para
frame_len = 128;
K = 256;
half_bin = K / 2 + 1;
D = frame_len;
win_len = K * 3;
comp_out = zeros(1, K);

% window didi
win_lpf_st = load('prototype_coef.mat');
win_lpf = win_lpf_st.prototype_filter_ori;
win = win_lpf';

% subband para
ana_win_echo = zeros(1, win_len);
ana_win_far = zeros(1, win_len);
sys_win = zeros(1, win_len);

% subband para
tap_low = 15;
tap_high = 10;
max_tap = max(tap_high, tap_low);
fre_low = 80;
fre_high = 81;
alpha_mse = zeros(half_bin,1);
alpha_mse(1:fre_low) = 0.90;
alpha_mse(fre_high:end) = 0.98;
alpha_mse_1 = 1-alpha_mse;

% background filter
subband_adf_num = zeros(1, half_bin);
subband_adf_num(1:fre_low) = tap_low;
subband_adf_num(fre_high:end) = tap_high;

subband_adf = zeros(half_bin, max_tap);
subband_adf_in = zeros(half_bin, max_tap);
subband_adf_mse = zeros(half_bin,1);

% foreground filter 
subband_fir = zeros(half_bin, max_tap);
subband_fir_in = zeros(half_bin, max_tap);
mic_in_mse = zeros(half_bin,1);
subband_fir_mse = zeros(half_bin,1);

% nlms para
myu = zeros(1, half_bin);
myu(1:fre_low) = 0.025;
myu(fre_high:end) = 0.025;
beta = 1e-4;

len_audio = min(length(y_echo),length(y_far));
tic;
for i = 1 : fix(len_audio/frame_len)
    % decompose echo subband
    in_frame_echo = y_echo((i - 1) * frame_len + 1 : i * frame_len);
    %ana_win_echo = [flipud(in_frame_echo)', ana_win_echo(1:end-frame_len)];
    ana_win_echo = [ana_win_echo(frame_len+1:end), in_frame_echo'];
    ana_win_echo_windowed = win' .* ana_win_echo;
    ana_wined_echo = ana_win_echo_windowed(1:K)+ana_win_echo_windowed(K+1:2*K)+ ana_win_echo_windowed(2*K+1:3*K);
    fft_out_echo = fft(ana_wined_echo, K);
    
    % decompose far subband
    in_frame_far = y_far((i - 1) * frame_len + 1 : i * frame_len);
    %ana_win_far = [flipud(in_frame_far)', ana_win_far(1:end-frame_len)];
    ana_win_far = [ana_win_far(frame_len + 1:end), in_frame_far'];
    ana_win_far_windowed = win' .* ana_win_far;
    ana_wined_far = ana_win_far_windowed(1:K)+ana_win_far_windowed(K+1:2*K)+ ana_win_far_windowed(2*K+1:3*K);
    fft_out_far = fft(ana_wined_far, K);
        
    % background filter out
    subband_adf_in(1:fre_low,1:tap_low) = [fft_out_far(1:fre_low)', subband_adf_in(1:fre_low,1:tap_low-1)];
    subband_adf_in(fre_high:end,1:tap_high) = [fft_out_far(fre_high:half_bin)', subband_adf_in(fre_high:end,1:tap_high-1)];
    subband_adf_out = sum(conj(subband_adf) .* subband_adf_in,2);
    subband_adf_err = fft_out_echo(1:half_bin)' - subband_adf_out;
    subband_adf_mse = alpha_mse .* subband_adf_mse + alpha_mse_1 .* subband_adf_err .* conj(subband_adf_err);

    % foreground filter out
    subband_fir_in(1:fre_low,1:tap_low) = [fft_out_far(1:fre_low)', subband_fir_in(1:fre_low,1:tap_low-1)];
    subband_fir_in(fre_high:end,1:tap_high) = [fft_out_far(fre_high:half_bin)', subband_fir_in(fre_high:end,1:tap_high-1)];    
    subband_fir_out = sum(conj(subband_fir) .* subband_fir_in,2);
    subband_fir_err = fft_out_echo(1:half_bin)' - subband_fir_out;
    subband_fir_mse = alpha_mse .* subband_fir_mse + alpha_mse_1 .* subband_fir_err .* conj(subband_fir_err);
    
    % mic mse
    mic_in_mse = alpha_mse .* mic_in_mse + alpha_mse_1 .* fft_out_echo(1:half_bin)' .* conj(fft_out_echo(1:half_bin))';
    
    % logic for update
    for j = 1: half_bin
        if subband_adf_mse(j) > mic_in_mse(j)*8
            subband_adf(j,:) = 0;
            subband_fir_mse(j) = 0;
            subband_adf_mse(j) = 0;
            mic_in_mse(j) = 0;
        elseif ((mic_in_mse(j) > subband_adf_mse(j) * 8) && (subband_adf_mse(j) < 0.5 * subband_fir_mse(j)))
            subband_fir(j,: ) = subband_adf(j,:);
            subband_fir_err(j) = subband_adf_err(j);
            subband_fir_mse(j) = 0;
            subband_adf_mse(j) = 0;
            mic_in_mse(j) = 0;
        end
        
        if subband_fir_mse(j) > mic_in_mse(j)*8
            subband_fir(j,:) = 0;
            subband_fir_mse(j) = 0;
            subband_adf_mse(j) = 0;
            mic_in_mse(j) = 0;
        elseif ((mic_in_mse(j) > subband_fir_mse(j) * 8) && (subband_fir_mse(j) < 0.5 * subband_adf_mse(j)))
            subband_adf(j,: ) = subband_fir(j,:);
            subband_adf_err(j) = subband_fir_err(j);    
            subband_fir_mse(j) = 0;
            subband_adf_mse(j) = 0;
            mic_in_mse(j) = 0;
        end
    end
    
    % nlms update
    norm = sum(subband_adf_in .* conj(subband_adf_in), 2) ./ subband_adf_num';
    alpha_nlms = myu' ./ (norm + beta);
    phi = alpha_nlms .* subband_adf_in .* conj(subband_adf_err);
    subband_adf = subband_adf + phi;
    
    ifft_in = [subband_fir_err', fliplr(conj(subband_fir_err(2:end-1)'))];
    
    % compose subband
    out = ifft(ifft_in, K) * K;
    %out = fliplr(out);
    win_in = [out, out, out];
    comp_out = win_in' .* win;
    sys_win = sys_win + comp_out';
    y_out((i-1) * frame_len + 1: i * frame_len) = sys_win(1 : frame_len) * frame_len ;
    sys_win = [sys_win(frame_len + 1 : end), zeros(1, frame_len)];

end
toc;
audiowrite(out_file, y_out/32678, fs_echo);
