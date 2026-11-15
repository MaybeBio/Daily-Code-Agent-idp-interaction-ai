#!/bin/sh
# UniPPI-OWL v1.8 Local Loopback Launcher (macOS / Linux / WSL)
# Binds strictly to 127.0.0.1:8766 for 100% private local browser execution.
PORT=8766
echo "======================================================================"
echo " UniPPI-OWL v1.8 Local Server running at: http://127.0.0.1:${PORT}/"
echo " Data Privacy Lock: 100% Local Loopback (Press Ctrl+C to stop)"
echo "======================================================================"
python3 -m http.server ${PORT} --bind 127.0.0.1
