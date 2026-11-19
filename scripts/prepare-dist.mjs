import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
const rootDir = path.resolve(__dirname, '..');
const distDir = path.join(rootDir, 'dist');

// Required frontend entrypoints
const REQUIRED_FILES = ['index.html'];

// Optional web data/assets in repository root
const OPTIONAL_FILES = [
  'mission_report.json',
  'favicon.ico',
  'favicon.svg',
  'favicon.png',
  'manifest.json',
];

// Optional static asset directories if added in the future
const OPTIONAL_DIRS = ['assets', 'public', 'static', 'css', 'js', 'img', 'images', 'fonts'];

function prepareDist() {
  fs.rmSync(distDir, { recursive: true, force: true });
  fs.mkdirSync(distDir, { recursive: true });

  const copied = [];

  for (const relPath of REQUIRED_FILES) {
    const src = path.join(rootDir, relPath);
    if (!fs.existsSync(src)) {
      throw new Error(`Required frontend file not found: ${relPath}`);
    }
    const dest = path.join(distDir, relPath);
    fs.copyFileSync(src, dest);
    copied.push(relPath);
  }

  for (const relPath of OPTIONAL_FILES) {
    const src = path.join(rootDir, relPath);
    if (fs.existsSync(src) && fs.statSync(src).isFile()) {
      const dest = path.join(distDir, relPath);
      fs.copyFileSync(src, dest);
      copied.push(relPath);
    }
  }

  for (const relDir of OPTIONAL_DIRS) {
    const src = path.join(rootDir, relDir);
    if (fs.existsSync(src) && fs.statSync(src).isDirectory()) {
      const dest = path.join(distDir, relDir);
      fs.cpSync(src, dest, { recursive: true });
      copied.push(`${relDir}/`);
    }
  }

  console.log(`[prepare:dist] Staged ${copied.length} item(s) into dist/: ${copied.join(', ')}`);
}

prepareDist();
