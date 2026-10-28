# node:20-slim (Debian glibc), NOT node:20-alpine. Next.js 16 + several
# of our deps (turbopack runtime, lucide, base-ui, plotly, molstar) ship
# native binaries whose musl builds either don't exist or have subtly
# different module shapes — the alpine container builds without error but
# fails at SSR runtime with "Cannot read properties of undefined (reading
# 'Provider')", because one of the named exports resolves to undefined
# only on musl. node:20-slim avoids the whole class of musl vs glibc bugs
# at the cost of ~100 MB extra image size.
FROM node:20-slim

WORKDIR /app

COPY package*.json ./
# `npm ci` (not `npm install`): refuses to mutate package-lock.json, so the
# container installs the EXACT same dep tree the user resolved on their Mac.
# Without this, `npm install` in the container can pick different versions
# for unpinned ranges and the resulting bundle mismatches what was tested.
RUN npm ci

COPY . .

# Next.js bakes NEXT_PUBLIC_* env vars at BUILD time — they cannot be
# changed at container-runtime. The frontend code reads
# `process.env.NEXT_PUBLIC_API_URL || "http://localhost:8000"`. Without
# this ENV the browser tries to fetch from the USER's localhost (not the
# server's). For prod we want same-origin requests so nginx routes /api/
# to the backend — pin the public domain here.
#
# Override at build with --build-arg NEXT_PUBLIC_API_URL=... for staging.
ARG NEXT_PUBLIC_API_URL=https://splisoforms.org
ENV NEXT_PUBLIC_API_URL=$NEXT_PUBLIC_API_URL

RUN npm run build

EXPOSE 3000

# Using Next config standalone is normally preferred, but this standard wrapper works for demonstration
CMD ["npm", "start"]
