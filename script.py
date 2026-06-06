import tkinter as tk
from tkinter import scrolledtext
from tkinter import simpledialog
try:
    import pyautogui
except Exception:
    pyautogui = None
import random
import time
import threading
import math


class HumanTypingApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Natural Code Typer")
        self.root.geometry("620x560")

        self.label = tk.Label(root, text="Paste your code below:", font=("Arial", 11))
        self.label.pack(pady=5)

        self.text_area = scrolledtext.ScrolledText(root, wrap=tk.WORD, width=70, height=20)
        self.text_area.pack(pady=10, padx=10, fill=tk.BOTH, expand=True)

        self.start_button = tk.Button(
            root, text="Start Typing (5s Delay)",
            command=self.start_typing_thread,
            bg="#4CAF50", fg="white", font=("Arial", 11, "bold")
        )
        self.start_button.pack(pady=10)

        self.status_label = tk.Label(root, text="Status: Ready", fg="gray")
        self.status_label.pack(pady=5)

        self.typing = False

    # ------------------------------------------------------------------ #
    #  Human-realism engine                                               #
    # ------------------------------------------------------------------ #

    def _human_delay(self, char, prev_char, base_delay, elapsed_budget, chars_remaining):
        """
        Multi-factor delay model that mimics real typing patterns.

        Factors:
          1. Base pacing – we keep a running "budget" so the total duration
             stays on target despite random variance.
          2. Digraph difficulty – some key-pairs are slow (e.g. 'th' fast,
             'qz' slow; same-finger rolls, long stretches).
          3. Post-punctuation micro-pause – humans pause after . , ; : ) }
          4. Word-boundary pause – slight hesitation at spaces.
          5. Line-boundary pause – bigger pause at newlines (thinking).
          6. Burst mode – occasional short runs of fast keys.
          7. Gaussian jitter – natural variance around the target delay.
          8. Rare long pause – simulates distraction / reading ahead.
        """

        # --- 1. Budget-aware base ---
        # Distribute remaining budget evenly, but allow ±40% jitter
        if chars_remaining > 0:
            budget_delay = elapsed_budget / chars_remaining
        else:
            budget_delay = base_delay
        target = max(0.01, budget_delay)

        # --- 2. Digraph difficulty ---
        difficulty = 1.0
        SLOW_PAIRS = {
            ('z','x'), ('x','z'), ('q','p'), ('p','q'),
            ('b','n'), ('n','b'), ('m','y'), ('y','m'),
        }
        FAST_PAIRS = {
            ('t','h'), ('h','e'), ('i','n'), ('e','r'),
            ('s','t'), ('i','o'), ('(',')')
        }
        if prev_char and (prev_char.lower(), char.lower()) in SLOW_PAIRS:
            difficulty *= random.uniform(1.4, 2.0)
        if prev_char and (prev_char.lower(), char.lower()) in FAST_PAIRS:
            difficulty *= random.uniform(0.5, 0.8)

        # --- 3. Post-punctuation pause ---
        if prev_char and prev_char in '.,:;!?':
            difficulty *= random.uniform(1.3, 1.9)

        # --- 4. Word-boundary pause ---
        if char == ' ':
            difficulty *= random.uniform(0.9, 1.4)

        # --- 5. Line-boundary (newline) pause ---
        if char in '\n\r':
            difficulty *= random.uniform(1.5, 3.0)

        # --- 6. Burst mode (5% chance) – run of fast keys ---
        if random.random() < 0.05:
            difficulty *= random.uniform(0.3, 0.6)

        # --- 7. Gaussian jitter (σ ≈ 20% of target) ---
        sigma = target * 0.20
        jitter = random.gauss(0, sigma)
        delay = max(0.008, target * difficulty + jitter)

        # --- 8. Rare long pause: 1.5% chance (reading / distraction) ---
        if random.random() < 0.015:
            delay += random.uniform(0.4, 1.8)

        return delay

    # ------------------------------------------------------------------ #
    #  Thread entry                                                        #
    # ------------------------------------------------------------------ #

    def start_typing_thread(self):
        if self.typing:
            return

        code = self.text_area.get("1.0", tk.END).strip()
        if not code:
            self._set_status("Status: No code provided!", 'red')
            return

        minutes = simpledialog.askfloat(
            "Typing Duration",
            "How many MINUTES should the typing last?",
            initialvalue=1.0, minvalue=0.1
        )
        if minutes is None:
            return

        self.target_duration = minutes * 60
        self.typing = True
        self._set_start_button(False)
        t = threading.Thread(target=self.simulate_typing, daemon=True)
        t.start()

    # ------------------------------------------------------------------ #
    #  Core typing loop                                                    #
    # ------------------------------------------------------------------ #

    def simulate_typing(self):
        code = self.text_area.get("1.0", tk.END).strip()
        if not code:
            self._set_status("Status: No code provided!", 'red')
            self._set_start_button(True)
            self.typing = False
            return

        total_chars = len(code)
        total_seconds = getattr(self, 'target_duration', 60)
        base_delay = total_seconds / total_chars if total_chars > 0 else 0.1

        info = (f"Target: {total_seconds:.1f}s | "
                f"Chars: {total_chars} | "
                f"Avg: {base_delay*1000:.0f}ms/char")

        for i in range(5, 0, -1):
            self._set_status(
                f"Status: Starting in {i}s... Switch to your editor!  {info}", 'blue')
            time.sleep(1)

        self._set_status("Status: Typing...", 'green')

        if pyautogui is None:
            self._set_status('pyautogui not available. Install it: pip install pyautogui', 'red')
            self._set_start_button(True)
            self.typing = False
            return

        try:
            pyautogui.PAUSE = 0
            pyautogui.FAILSAFE = True   # slam mouse to any corner to abort
        except Exception:
            pass

        start_time = time.time()
        prev_char = None

        for i, char in enumerate(code):
            chars_remaining = total_chars - i
            elapsed = time.time() - start_time
            time_remaining = total_seconds - elapsed
            budget = time_remaining   # seconds left to spend

            delay = self._human_delay(
                char, prev_char, base_delay,
                elapsed_budget=budget,
                chars_remaining=chars_remaining
            )

            # Sleep BEFORE typing (reaction-time model)
            if delay > 0:
                time.sleep(delay)

            try:
                if char in ('\n', '\r'):
                    pyautogui.press('enter')
                elif char == '\t':
                    pyautogui.press('tab')
                else:
                    pyautogui.write(char, interval=0)
            except Exception as e:
                self._set_status(f"Typing aborted: {e}", 'red')
                break

            prev_char = char

            # Live progress every ~50 chars
            if i % 50 == 0:
                pct = int(100 * i / total_chars)
                elapsed_now = time.time() - start_time
                self._set_status(
                    f"Status: Typing… {pct}% | {elapsed_now:.0f}s elapsed", 'green')

        elapsed = time.time() - start_time
        self._set_status(f"Status: Done! Took {elapsed:.1f}s of {total_seconds:.1f}s target.", 'gray')
        self._set_start_button(True)
        self.typing = False

    # ------------------------------------------------------------------ #
    #  Helpers                                                             #
    # ------------------------------------------------------------------ #

    def _set_status(self, text, fg='gray'):
        try:
            self.root.after(0, lambda: self.status_label.config(text=text, fg=fg))
        except Exception:
            pass

    def _set_start_button(self, enabled):
        state = 'normal' if enabled else 'disabled'
        try:
            self.root.after(0, lambda: self.start_button.config(state=state))
        except Exception:
            pass


if __name__ == "__main__":
    root = tk.Tk()
    app = HumanTypingApp(root)
    root.mainloop()