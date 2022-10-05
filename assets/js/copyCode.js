/*Adapted from https://www.aleksandrhovhannisyan.com/blog/how-to-add-a-copy-to-clipboard-button-to-your-jekyll-blog/*/
const copyCodeButtons = document.querySelectorAll('.copy-code-button');
copyCodeButtons.forEach(copyCodeButton) => {
  const code = copyCodeButton.previousSibling.innerText;

  copyCodeButton.addEventListener('click', () => {
    window.navigator.clipboard.writeText(code);
    /*window.alert(code);*/
    copyCodeButton.classList.remove('fa-copy');
    copyCodeButton.classList.add('fa-check');

    setTimeout(() => {
      copyCodeButton.classList.remove('fa-check');
      copyCodeButton.classList.add('fa-copy');
    }, 2000);
  });
});

/* ORIGINAL CODE, FOR REFERENCE
const codeBlocks = document.querySelectorAll('.code-header + .highlighter-rouge');
const copyCodeButtons = document.querySelectorAll('.copy-code-button');

copyCodeButtons.forEach((copyCodeButton, index) => {
  const code = codeBlocks[index].innerText;

  copyCodeButton.addEventListener('click', () => {
    // Copy the code to the user's clipboard
    window.navigator.clipboard.writeText(code);

    // Update the button text visually
    const { innerText: originalText } = copyCodeButton;
    copyCodeButton.innerText = 'Copied!';

    // (Optional) Toggle a class for styling the button
    copyCodeButton.classList.add('copied');

    // After 2 seconds, reset the button to its initial UI
    setTimeout(() => {
      copyCodeButton.innerText = originalText;
      copyCodeButton.classList.remove('copied');
    }, 2000);
  });
});*/
