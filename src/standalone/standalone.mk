include src/chemistry/mk/core.mk

# Add these sources only to a standalone host build. The chemistry
# core source list above must not contain host adapters.
CHEMISTRY_STANDALONE_ADAPTER_SOURCES += \
  src/standalone/chemistry_standalone_adapter.c \
  src/standalone/chemistry_standalone_background.c
