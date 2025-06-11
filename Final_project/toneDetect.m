function [maxNote, tone, status, status_fft, tone_fft] = toneDetect(xx, fs, DEBUG, filename)
    if (nargin < 3)
        DEBUG = 0;
    end
    f0 = 261.63;
    fprintf("fs = %d\n",fs);
    %fs = 8000;
    all = xx';
    maxpower = 0;
    maxNote = -1;
    keys = 0:12;
    values = {
        "C4", "C#4", "D4", "D#4", "E4", "F4", "F#4", "G4", "G#4", "A4", "A#4", "B4", "C5"
    };
    notemap = containers.Map(keys, values);
    freq_map = [ ...
        261.63, 277.18, 293.66, 311.13, 329.63, ...
        349.23, 369.99, 392.00, 415.30, 440.00, ...
        466.16, 493.88, 523.25];
    powerlist = zeros(1, 13);
    bw = 10;
    for note = 0:12
        freq = f0 * 2^(note / 12);
        [b, a] = iirpeak(freq / (fs/2), bw / (fs/2));
        yy = filter(b, a, xx);
        normFactor = sum(filter(b, a, ones(size(xx))).^2);
        power = sum(yy.^2) / normFactor;
        powerlist(note + 1) = power;
        if power > maxpower
            maxpower = power;
            maxNote = note;
        end
        all = [all yy'];
    end

    tone = notemap(maxNote);

    N = length(xx);
    X = abs(fft(xx));
    f_axis = (0:N-1)*(fs/N);
    half_X = X(1:floor(N/2));
    half_f = f_axis(1:floor(N/2));
    [~, idx] = max(half_X(:));
    f_peak = half_f(idx);
    [~, closest_idx] = min(abs(freq_map - f_peak));
    closest_freq = freq_map(closest_idx);
    tone_fft = values{closest_idx};
    diff_Hz = abs(f_peak - closest_freq);

    if f_peak < freq_map(1) || f_peak > freq_map(end)
        status_fft = "Out of Range";
    elseif diff_Hz <= 5
        status_fft = "Exact Match";
    elseif diff_Hz <= 15
        status_fft = "Detuned";
    else
        status_fft = "Out of Range";
    end

    normalized = powerlist / sum(powerlist + eps);
    entropy = -sum(normalized .* log(normalized + eps));
    %normalize_maxpower = max(powerlist) / sum(powerlist + eps);
    sorted_normalized_power = sort(normalized, 'descend');
    top1_energy_ratio = sorted_normalized_power(1);

    if top1_energy_ratio > 0.7 && entropy < 1.0
        status = "Exact Match";
    elseif top1_energy_ratio < 0.4 || entropy > 1.8
        status = "Out of Range";
    else
        status = "Detuned";
    end

    if DEBUG
        figure;
        bar(normalized, 'FaceColor', [0.2 0.6 0.8]);
        title(sprintf('Note Energy Distribution — Detected: %s (%s)', tone, status), 'FontSize', 12);
        xlabel('Note');
        ylabel('Power (Energy)');
        xticks(1:13);
        xticklabels(values);
        xtickangle(45);
        grid on;

        fprintf("%s | FFT: %s | Note: %s | Reason: peak at %.2f Hz\n", ...
            filename, status_fft, tone_fft, f_peak);

        fprintf("%s | Power: %s | Note: %s | Reason: max energy at filter %s\n", ...
            filename, status, tone, tone);
        %soundsc(all, fs);
    end
end
