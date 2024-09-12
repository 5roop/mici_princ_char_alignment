from os import environ

try:
    wavpath = snakemake.input.wav
    jsonpath = snakemake.input.json
    outputpath = snakemake.output.json
    environ["CUDA_VISIBLE_DEVICES"] = snakemake.params.cuda
except NameError:
    wavpath = "data/MPwav/MP_00.wav"
    jsonpath = "data/MPjson/MP_00.json"
    outputpath = "brisi.json"
    environ["CUDA_VISIBLE_DEVICES"] = "4"

from pathlib import Path
import json
import numpy as np
from transformers import Wav2Vec2Processor, AutoModelForCTC
from datasets import Dataset, Audio
import torch

w2v2_model = "classla/wav2vec2-large-slavic-parlaspeech-hr"
processor = Wav2Vec2Processor.from_pretrained(str(w2v2_model))
model = AutoModelForCTC.from_pretrained(str(w2v2_model)).cuda()


sampling_rate = 16000
audio = Dataset.from_dict({"audio": [wavpath]}).cast_column(
    "audio", Audio(sampling_rate=sampling_rate, mono=True)
)[0]["audio"]["array"]
mp = json.loads(Path(jsonpath).read_text())


def get_char_alignment(trans, audio):
    labels = [
        processor.tokenizer.encoder.get(i, processor.tokenizer.unk_token_id)
        for i in trans.casefold()
    ]
    processed = processor.feature_extractor(
        audio, sampling_rate=sampling_rate, return_tensors="pt"
    )

    with torch.no_grad():
        logits = (
            model(
                input_values=processed["input_values"].cuda(),
                attention_mask=processed["attention_mask"].cuda(),
            )
            .logits.cpu()
            .numpy()
            .squeeze()
        )

    def viterbi(
        logits: np.ndarray,
        labels: np.ndarray,
        pad_id: int = processor.tokenizer.pad_token_id,
    ) -> tuple[list[int], np.ndarray]:
        T = logits.shape[0]
        N = len(labels)
        delta = np.zeros((T, N))
        psi = np.zeros((T, N), dtype=int)

        delta[0, :] = -np.inf
        delta[0, 0] = logits[0][pad_id]

        a = np.zeros((N,))
        b = np.zeros((N,))
        o = np.zeros((N,))
        M = np.arange(N)

        for t in range(1, T):
            delta[t, :] = -np.inf
            beg = np.clip(t - T + N, 0, N - 1)
            end = np.clip(t + 1, 0, N)

            r = np.arange(beg, end)
            ra = np.clip(r, 1, N - 1)

            a[ra] = delta[t - 1][ra - 1]
            b[r] = delta[t - 1][r]

            o[r] = logits[t][labels[r]]
            opad = logits[t][pad_id]

            m = np.zeros(N, dtype=bool)
            m[ra] = (a > b)[ra]
            delta[t][m] = a[m] + o[m]
            psi[t][m] = M[m] - 1

            m = np.zeros(N, dtype=bool)
            m[r] = (a <= b)[r]
            delta[t][m] = b[m] + opad
            psi[t][m] = M[m]

        bestpath = [delta[T - 1].argmax()]
        for t in range(T - 2, -1, -1):
            bestpath.append(psi[t + 1][bestpath[-1]])
        bestpath.reverse()
        return bestpath, delta

    def align(
        logits: np.ndarray, labels: list[int], fr_len: float = 0.02
    ) -> np.ndarray:
        # Use character labels directly (no word-level ID handling)
        bestpath, _ = viterbi(logits, np.array(labels))

        # Calculate the duration for each aligned character
        char_durations = []
        current_char = bestpath[0]  # Start with the first character
        start_frame = 0

        for t, p in enumerate(bestpath):
            if p != current_char:  # If we encounter a new character
                char_durations.append(
                    (start_frame, t)
                )  # Append start and end frame for the last character
                current_char = p  # Update to the new character
                start_frame = t  # Reset the start frame

        # Append the duration for the final character
        char_durations.append((start_frame, len(bestpath)))

        # Each duration corresponds to a character in 'trans' and should have the same length
        assert len(char_durations) == len(
            labels
        ), f"Mismatch: {len(char_durations)} durations, {len(labels)} labels."

        # Convert frame durations into time durations using fr_len (frame length)
        return np.array([[start, end] for start, end in char_durations]) * fr_len

    # Align characters and get the result
    char_lens = align(logits, labels)
    return char_lens


results = dict()
for i, entry in enumerate(mp):
    start = entry["time_s"]
    end = entry["time_e"]
    start_index = int(start * sampling_rate)
    end_index = int(end * sampling_rate)
    trans = entry["text"]
    audio_segment = audio[start_index:end_index]
    char_lens = get_char_alignment(trans, audio_segment)
    char_lens = (char_lens + start).round(2).tolist()
    word_lens = [[i["time_s"], i["time_e"]] for i in entry["words"]]
    for j, crow in enumerate(char_lens):
        s, e = crow
        for ws, we in word_lens:
            if (s > ws) and (e > we) and (s < we) and (e > ws):
                char_lens[j] = [s, we]
                break
    first_word_offset = entry["words"][0]["char_s"]
    current_results = []
    for j, (time_s, time_e) in enumerate(char_lens):
        assert time_s <= time_e, "Negative duration!"
        current_results.append(
            {
                "char": trans[j],
                "char_s": first_word_offset + j,
                "char_e": first_word_offset + j + 1,
                "time_s": time_s,
                "time_e": time_e,
            }
        )
    assert len(current_results) == len(trans)
    results[i] = current_results
for i, _ in enumerate(mp):
    mp[i]["chars"] = results[i]
Path(outputpath).write_text(json.dumps(mp, ensure_ascii=False, indent=4))
