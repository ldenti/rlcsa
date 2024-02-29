#include <cstdint>
#include <getopt.h>
#include <iostream>
#include <string>
#include <zlib.h>

#include "kseq.h"
#include "rlcsa.h"
#include "rlcsa_builder.h"

KSEQ_INIT(gzFile, gzread)

static unsigned char seq_nt6_table[128] = {
    0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1,
    5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};

static inline uint kputsn(const char *p, uint l, kstring_t *s) {
  if (s->l + l + 1 >= s->m) {
    char *tmp;
    s->m = s->l + l + 2;
    kroundup32(s->m);
    if ((tmp = (char *)realloc(s->s, s->m)))
      s->s = tmp;
    else
      return EOF;
  }
  memcpy(s->s + s->l, p, l);
  s->l += l;
  s->s[s->l] = 0;
  return l;
}

int main(int argc, char **argv) {
  double rt = realtime();
  uint sample_rate = 0; // (1 << 31);
  int block_size = 32;
  bool reverse = false;
  int threads = 1;
  std::string index_prefix = "RLCSA";

  int c;
  while ((c = getopt(argc, argv, "i:@:rvh")) >= 0) {
    switch (c) {
    case 'i':
      index_prefix = optarg;
      continue;
    case 'r':
      reverse = true;
      continue;
    case '@':
      threads = atoi(optarg);
      continue;
    default:
      return 1;
    }
  }
  if (argc - optind < 1) {
    return 1;
  }

  CSA::RLCSABuilder builder(block_size, sample_rate, 0, threads, NULL);

  int64_t m = (int64_t)(.97 * 10 * 1024 * 1024 * 1024) + 1;
  gzFile fp;
  kseq_t *ks;
  uint8_t *s;
  kstring_t buf = {0, 0, 0};
  int i, l;
  char *fa_path;
  while (optind < argc) {
    fa_path = argv[optind++];

    fp = gzopen(fa_path, "rb");
    ks = kseq_init(fp);
    while ((l = kseq_read(ks)) >= 0) {
      s = (uint8_t *)ks->seq.s;

      // change encoding
      for (i = 0; i < l; ++i)
        s[i] = s[i] < 128 ? seq_nt6_table[s[i]] : 5;

      // Add forward to buffer
      kputsn((char *)ks->seq.s, ks->seq.l + 1, &buf);
      if (reverse) {
        // Add reverse to buffer
        for (i = 0; i < (l >> 1); ++i) {
          int tmp = s[l - 1 - i];
          tmp = (tmp >= 1 && tmp <= 4) ? 5 - tmp : tmp;
          s[l - 1 - i] = (s[i] >= 1 && s[i] <= 4) ? 5 - s[i] : s[i];
          s[i] = tmp;
        }
        if (l & 1)
          s[i] = (s[i] >= 1 && s[i] <= 4) ? 5 - s[i] : s[i];
        kputsn((char *)s, ks->seq.l + 1, &buf);
      }
      if (buf.l >= m) {
        CSA::RLCSA *index =
            new CSA::RLCSA((CSA::uchar *)buf.s, buf.l, block_size, sample_rate,
                           threads, false);
        fprintf(stderr,
                "[M::%s] created partial index from %ld symbols in %.3f sec\n",
                __func__, (long)buf.l, realtime() - rt);
        rt = realtime();

        index->writeTo(fa_path); // TODO: store to tmp dir
        delete index;
        fprintf(stderr, "[M::%s] stored partial index in %.3f sec\n", __func__,
                realtime() - rt);
        rt = realtime();
        builder.insertFromFile(fa_path, (CSA::uchar *)buf.s);
        fprintf(stderr, "[M::%s] merged partial index in %.3f sec\n", __func__,
                realtime() - rt);
        rt = realtime();

        buf.l = 0;
      }
    }
    if (buf.l) {
      CSA::RLCSA *index = new CSA::RLCSA((CSA::uchar *)buf.s, buf.l, block_size,
                                         sample_rate, threads, false);
      fprintf(stderr,
              "[M::%s] created partial index from %ld symbols in %.3f sec\n",
              __func__, (long)buf.l, realtime() - rt);

      index->writeTo(fa_path); // TODO: store to tmp dir
      delete index;
      fprintf(stderr, "[M::%s] stored partial index in %.3f sec\n", __func__,
              realtime() - rt);
      builder.insertFromFile(fa_path, (CSA::uchar *)buf.s);
      fprintf(stderr, "[M::%s] merged partial index in %.3f sec\n", __func__,
              realtime() - rt);

      buf.l = 0;
    }
    kseq_destroy(ks);
    gzclose(fp);
  }
  free(buf.s);
  CSA::RLCSA *rlcsa = builder.getRLCSA();
  if (rlcsa != 0 && rlcsa->isOk()) {
    // std::cout << std::endl;
    // rlcsa->printInfo();
    // rlcsa->reportSize(true);
    rlcsa->writeTo(index_prefix);
  }
  delete rlcsa;
  fprintf(stderr,
                "[M::%s] Stored full index in %.3f sec\n",
                __func__, realtime() - rt);
  return 0;
}
