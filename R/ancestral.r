# mispolarized (HC)
sites_mp_hc <- full_data$sites %>%
  filter(ALT==AA)

# mispolarized (HC+LC)
sites_mp_all <- full_data$sites %>%
  filter(tolower(ALT)==tolower(AA))

# correct (HC)
sites_c_hc <- full_data$sites %>%
  filter(REF==AA)

# correct (HC+LC)
sites_c_all <- full_data$sites %>%
  filter(tolower(REF)==tolower(AA))

#
sites_c_lc <- full_data$sites %>%
  filter(tolower(REF)==AA | AA %in% c("-", "N", "."))


  # grepl("^[[:upper:]]+$", s)
