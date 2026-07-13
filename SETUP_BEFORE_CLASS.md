# Setup — do this BEFORE the course

Please do these 4 steps **at home, on a good internet connection**, before Day 1. The
downloads are several GB — you do not want to run them on classroom Wi-Fi. Budget ~30–60
minutes (mostly unattended downloading).

If anything fails, message us with the error text and your OS — we'll sort it out before
class rather than during it.

---

## 1. Install Docker Desktop

- Download & install from **https://docs.docker.com/desktop**
  (macOS Intel/Apple-Silicon, Windows 10/11, or Linux).
- On **Windows**, enable **WSL2** when prompted (Docker's installer guides you).
- **Start Docker Desktop** and leave it running. Check it works:

```bash
docker run --rm hello-world
```
You should see "Hello from Docker!".

> Give Docker enough resources: Docker Desktop → **Settings → Resources** →
> at least **6 GB RAM** and **4 CPUs**, and make sure you have **~20 GB free disk**.

---

## 2. Get the course code

```bash
git clone https://github.com/DomeJoyce/epicrops-course-2026.git
cd epicrops-course-2026
```
(No git? Download the ZIP from the GitHub page and unzip it.)

---

## 3. Get the Day 1 image (lncRNA/mRNA)

This one image contains all Day-1 tools **and** the data (~6 GB on disk):

```bash
docker pull leogiuffre/lncrna-mnps-workshop:1.0
```
> **Apple-Silicon Macs (M1/M2/M3):** `docker pull --platform linux/amd64 leogiuffre/lncrna-mnps-workshop:1.0`

Quick check it starts (then press **Ctrl-C** to stop it):
```bash
docker compose up day1
```
When you see a JupyterLab URL / `http://localhost:8888`, it works — stop it with Ctrl-C.

---

## 4. Get the Day 2 image (WGBS)

**Pull the pre-built image** (fast — a few minutes):
```bash
docker pull djoyce86/epi-code-practical:2026
docker tag djoyce86/epi-code-practical:2026 epi-code-practical:latest
```
Then check the tools are present:
```bash
docker compose run --rm course bash -lc 'bash $SCRIPTS_DIR/validate_env.sh'
```
It should finish with **`ALL GOOD — environment is ready.`**

> *No internet at the venue, or prefer building it yourself? It also works from source
> (~20–40 min the first time):* `docker compose build course`

---

## You're ready ✅

Both images pulled/built and both checks green → you're set for the course. The
WGBS **read download** for Day 2 (`download_data.sh`) is quick and we'll do it together in
class, so you don't need to run it now.

See you in Florence!
